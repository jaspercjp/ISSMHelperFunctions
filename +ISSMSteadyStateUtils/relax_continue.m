function [md, provenance] = relax_continue(config_path)
%RELAX_CONTINUE  Evolve an existing final_md.mat another N years, on the SAME
%   mesh, with the SAME physics, writing back IN PLACE.
%
%   [MD, PROVENANCE] = ISSMSteadyStateUtils.relax_continue(CONFIG_PATH)
%
%   One argument, exactly like calc_steady_state. CONFIG_PATH is the derived
%   config that Scripts/Workflows/SteadyState/relax.py materialises under
%   $SCRATCH/relax_configs/ -- the tracked repo configs are never edited.
%
% ===========================================================================
% THE PROBLEM THIS SOLVES
% ===========================================================================
% calc_steady_state.m runs a fixed four-phase pipeline and final_md.mat spans
% only AMR_INTERP_SETTINGS.NUM_YRS/2 -- currently 250 yr for flatbed/constslope.
% That is far short of convergence: the flatbed seed backing the shipped figures
% relaxed for ~1700 yr in manually assembled continuation chunks, and x_gl was
% STILL migrating at ~0.005 km/yr there. Before this function the only way to get
% more relax was to raise NUM_YRS and re-run from a cold bamg mesh -- throwing
% away hours of correct compute and, worse, producing a DIFFERENT mesh and
% therefore a different model, which every downstream IceMelt experiment would
% silently inherit.
%
% ===========================================================================
% WHAT MAKES LEG N+1 THE SAME EXPERIMENT  (design §5)
% ===========================================================================
% Carried over untouched, simply by loading final_md.mat and not touching it:
%   md.friction.*            (law class and coefficients)
%   md.materials.*           (rheology B, n)
%   md.stressbalance.spc*    (BCs, including the NoSlip walls)
%   md.smb.mass_balance, md.mesh.*, md.geometry.bed,
%   md.masstransport.stabilization, md.levelset.stabilization,
%   md.groundingline.*, md.transient.is*,
%   md.transient.requested_outputs, md.settings.solver_residue_threshold,
%   md.flowequation (SSA)
%
% Set by transientrestart, from the seed snapshot:
%   md.timestepping.start_time, md.initialization.{vx,vy,vel,pressure},
%   md.geometry.{base,thickness,surface}, md.mask.ocean_levelset
%
% Set HERE, and only these:
%   md.timestepping.final_time, md.miscellaneous.name, md.cluster,
%   md.miscellaneous.dummy.relax_provenance
%
% DELIBERATELY NOT CALLED -- each appears in calc_steady_state's relax phases and
% each would break the continuation:
%   ISSMModelUtils.ReinitializeAMRModel  rebuilds md from the domain outline and
%                                        discards the refined mesh
%   PARAM_FUNC (ParamMISMIP/ParamConstSlope)
%                                        re-derives geometry/friction from the
%                                        config, overwriting the evolved state
%   ISSMParamUtils.init_BC               re-imposes BCs on a mesh it did not build
%   ISSMMeshUtils.copy_amr_mesh          only meaningful after a remesh; none here
%   bamg / any remesh                    violates the mesh-verbatim contract that
%                                        every IceMelt experiment depends on
%
% BIT-REPRODUCIBILITY, HONESTLY. relax_continue twice is NOT bit-identical to
% calc_steady_state with a 500 yr final relax, for two inherent reasons:
% transientrestart reinitialises the velocity solve from the STORED Vx/Vy rather
% than the solver's in-memory iterate, and the leg boundary re-anchors the
% timestep sequence on a snapshot time. The difference is at the level of one
% nonlinear-solve convergence tolerance (solver_residue_threshold, 1e-6 here),
% i.e. far below the MS_diff_H ~ 4e-2 signal being measured. "Bit-defensible"
% here means "provably the same model and the same mesh, restarted at a recorded
% time" -- NOT "bit-identical trajectory". Do not claim more than that.
%
% ===========================================================================
% ORDER OF OPERATIONS
% ===========================================================================
%    1  read config; resolve RUN_DIR; getenv DATA_ROOT / ISSM_EXECUTION_ROOT
%    2  acquire RUN_DIR/relax.lock                                    -- G4
%    3  refuse if RUN_DIR/final_md.tmp-*.mat exists                   -- G3
%    4  load final_md.mat
%    5  PRE-FLIGHT ASSERTIONS                          -- G1 G2 G5 G7 G8 G11
%    6  read existing provenance from md.miscellaneous.dummy
%    7  keep a handle on the existing history
%    8  md = transientrestart(md)
%    9  final_time = start_time + ADD_YEARS
%   10  name / cluster
%   11  solve
%   12  POST-SOLVE ASSERTIONS                              -- G10 G12 G13
%   13  merge + thin                                            -- G6
%   14  provenance
%   15  save_model_atomic                                       -- G14
%   16  sidecar JSON; release lock
%
% Every guard is FATAL BEFORE THE SOLVE unless marked otherwise. An 8-hour solve
% that then refuses to write is the worst possible outcome, so the checks that
% can be made early are made early.

    t_wall_start = tic;

    % ------------------------------------------------------------- 1. config --
    if ~exist(config_path, 'file')
        error('ISSMSteadyStateUtils.relax_continue: no such config: %s', config_path);
    end
    config = jsondecode(fileread(config_path));

    data_root = getenv('DATA_ROOT');
    if isempty(data_root)
        error('ISSMHelperFunctions:noDataRoot', 'DATA_ROOT is not set; the launchers (Core.issm_run, _calc_steady_state.sh) export it.');
    end
    exec_root = getenv('ISSM_EXECUTION_ROOT');
    if isempty(exec_root)
        error('ISSMHelperFunctions:noExecRoot', 'ISSM_EXECUTION_ROOT is not set; the launchers (Core.issm_run, _calc_steady_state.sh) export it.');
    end

    run_dir = char(config.RUN_DIR);
    if run_dir(end) == filesep; run_dir(end) = []; end
    add_years = double(config.ADD_YEARS);
    keep_generations = 1;
    if isfield(config, 'KEEP_GENERATIONS') && ~isempty(config.KEEP_GENERATIONS)
        keep_generations = double(config.KEEP_GENERATIONS);
    end
    num_proc = 4;
    if isfield(config, 'NUM_PROC') && ~isempty(config.NUM_PROC)
        num_proc = double(config.NUM_PROC);
    end
    if ~isfield(config, 'RETENTION')
        error(['ISSMSteadyStateUtils.relax_continue: the config has no ' ...
               'RETENTION block. The policy''s source of truth is ' ...
               'Core.relax_retention.STANDARD_POLICY and relax.py bakes it in; ' ...
               'a config without it was hand-written.']);
    end
    policy = config.RETENTION;

    model_path = fullfile(run_dir, 'final_md.mat');

    fprintf('==============================================================\n');
    fprintf(' relax_continue\n');
    fprintf('   config     : %s\n', config_path);
    fprintf('   run_dir    : %s\n', run_dir);
    fprintf('   add_years  : %g\n', add_years);
    fprintf('   retention  : %s\n', policy.profile);
    fprintf('   generations: %d\n', keep_generations);
    fprintf('   host/job   : %s / %s\n', oshostname(), getenv('SLURM_JOB_ID'));
    fprintf('==============================================================\n');

    % ------------------------------------------------------- G11 (path sanity) --
    % A typo'd RUN_DIR must not be able to destroy an unrelated file.
    expected_parent = fullfile(data_root, 'steady_states');
    if ~startsWith(run_dir, expected_parent)
        error(['ISSMSteadyStateUtils.relax_continue [G11]: RUN_DIR\n  %s\nis ' ...
               'not under\n  %s\nThis function REWRITES final_md.mat in place; ' ...
               'it will not do that outside the steady-state tree.'], ...
              run_dir, expected_parent);
    end
    if ~exist(model_path, 'file')
        error(['ISSMSteadyStateUtils.relax_continue [G11]: %s does not exist. ' ...
               'A continuation needs something to continue.'], model_path);
    end

    % ---------------------------------------------------------- G3 (stale tmp) --
    tmps = dir(fullfile(run_dir, 'final_md.tmp-*.mat'));
    if ~isempty(tmps) && ~(isfield(config, 'FORCE') && config.FORCE)
        names = strjoin({tmps.name}, ', ');
        error(['ISSMSteadyStateUtils.relax_continue [G3]: %s contains ' ...
               'orphan write temporaries: %s\nA previous leg died mid-write. ' ...
               'final_md.mat itself is intact (save_model_atomic never touches ' ...
               'it until the tmp verifies), but SOMEONE MUST LOOK before it is ' ...
               'overwritten again. Delete the tmp file(s) by hand, or pass ' ...
               '--force.'], run_dir, names);
    end

    % ------------------------------------------------------------- 2. G4 lock --
    lock_path = fullfile(run_dir, 'relax.lock');
    acquire_lock(lock_path, config);
    cleanup = onCleanup(@() release_lock(lock_path));

    % --------------------------------------------------------------- 4. load --
    fprintf('[relax_continue] loading %s ...\n', model_path);
    loaded = load(model_path, 'md');
    md = loaded.md;
    clear loaded

    nvert_in    = md.mesh.numberofvertices;
    mesh_x_in   = md.mesh.x;
    smb_mean_in = mean(md.smb.mass_balance);
    dt_in       = md.timestepping.time_step;

    % ------------------------------------------------- 5. pre-flight assertions --

    % G1 -- a continuation that remeshes would silently change the resolution of
    % every downstream melt run. That is the exact failure the uniform 1 km
    % campaign exists to prevent, so it is asserted, never assumed.
    if md.transient.amr_frequency ~= 0
        error(['ISSMSteadyStateUtils.relax_continue [G1]: ' ...
               'md.transient.amr_frequency = %g, expected 0. final_md.mat is ' ...
               'written by the MESH-FIXED relax. If a continuation ran with AMR ' ...
               'on, the mesh would change and every IceMelt experiment seeded ' ...
               'off this run would silently inherit a different resolution.'], ...
              md.transient.amr_frequency);
    end

    % G5 -- a stray TransientSolution2 means a previous transientrestart was never
    % cleaned up. Merging on top of that would double-count the history.
    res_fields = fieldnames(md.results);
    if ~(numel(res_fields) == 1 && strcmp(res_fields{1}, 'TransientSolution'))
        error(['ISSMSteadyStateUtils.relax_continue [G5]: md.results holds {%s}, ' ...
               'expected exactly TransientSolution. A leftover ' ...
               'TransientSolution2 is what transientrestart leaves behind when ' ...
               'ISSMMiscUtils.delete_prev_transientsol was not run; merging on ' ...
               'top of it would double-count the whole history.'], ...
              strjoin(res_fields', ', '));
    end

    history = md.results.TransientSolution;
    n_hist  = numel(history);
    if n_hist < 1
        error('ISSMSteadyStateUtils.relax_continue: TransientSolution is empty.');
    end

    % G2 -- the seed record must be COMPLETE. transientrestart reads
    % TransientSolution(end) and nothing else; if someone thinned too hard, it
    % would install empty velocity/geometry fields and the leg would start from
    % garbage -- and it would do so SILENTLY, because isfield() is true on a
    % blanked field (see thin_transientsol.m, "THE MATLAB CONSTRAINT").
    seed_required = {'time', 'Vx', 'Vy', 'Vel', 'Pressure', 'Base', ...
                     'Thickness', 'MaskOceanLevelset'};
    seed = history(end);
    missing = {};
    for k = 1:numel(seed_required)
        f = seed_required{k};
        if ~isfield(seed, f) || isempty(seed.(f))
            missing{end+1} = f; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error(['ISSMSteadyStateUtils.relax_continue [G2]: the seed record ' ...
               'TransientSolution(end) is missing or EMPTY in: %s. ' ...
               'transientrestart reads only that record, so the leg would ' ...
               'start from empty initialisation fields. This file has been ' ...
               'over-thinned and cannot be continued.'], strjoin(missing, ', '));
    end

    % G7 -- the `--years 0` typo, and a leg too short to produce snapshots.
    if ~(add_years > 0)
        error(['ISSMSteadyStateUtils.relax_continue [G7]: ADD_YEARS = %g. ' ...
               'A continuation of zero years is a typo.'], add_years);
    end
    if add_years < 10 * dt_in
        error(['ISSMSteadyStateUtils.relax_continue [G7]: ADD_YEARS = %g is ' ...
               'under 10 timesteps at dt = %g. Too short to produce a usable ' ...
               'snapshot series.'], add_years, dt_in);
    end

    % Optional, LOUD overrides. Changing dt between legs makes the merged series
    % a splice of two numerical schemes; it is not bit-defensible and this
    % function will not pretend otherwise. Recorded per-leg in provenance.
    if isfield(config, 'TIME_STEP') && ~isempty(config.TIME_STEP)
        fprintf(2, ['[relax_continue] *** TIME_STEP OVERRIDE %g -> %g. The ' ...
                    'merged series now spans TWO distinct time steps and is a ' ...
                    'splice of two numerical schemes. Recorded in ' ...
                    'provenance.legs[].time_step. ***\n'], ...
                dt_in, config.TIME_STEP);
        md.timestepping.time_step = config.TIME_STEP;
    end
    if isfield(config, 'OUTPUT_FREQ') && ~isempty(config.OUTPUT_FREQ)
        fprintf('[relax_continue] output_frequency override %g -> %g\n', ...
                md.settings.output_frequency, config.OUTPUT_FREQ);
        md.settings.output_frequency = config.OUTPUT_FREQ;
    end

    % G8 -- provenance integrity. NOT fatal: the whole campaign is currently
    % running against an uncommitted submodule working tree. Recorded as
    % helpers_dirty so a run made against unpushed code is identifiable later.
    [helpers_sha, helpers_dirty] = helpers_git_state();
    if helpers_dirty
        fprintf(2, ['[relax_continue] WARNING [G8]: the ISSMHelperFunctions ' ...
                    'working tree is DIRTY at %s. This leg is not reproducible ' ...
                    'by anyone else. Recorded as helpers_dirty = true.\n'], ...
                helpers_sha);
    end

    % ------------------------------------------------------ 6. read provenance --
    prov = existing_provenance(md, run_dir);

    fprintf('[relax_continue] history: %d snapshots, t = %.4f -> %.4f\n', ...
            n_hist, double(history(1).time), double(history(end).time));

    % ------------------------------------------------------------ 8-10. restart --
    md = transientrestart(md);
    % transientrestart moves the array into TransientSolution2 and empties
    % TransientSolution. We already hold `history`, so drop its copy immediately
    % -- leaving it would double the resident footprint AND trip G5 next leg.
    if isfield(md.results, 'TransientSolution2')
        md.results = rmfield(md.results, 'TransientSolution2');
    end

    start_time = md.timestepping.start_time;
    md.timestepping.final_time = start_time + add_years;
    leg_number = numel(prov.legs) + 1;

    date_str = ISSMMiscUtils.get_date_str();
    % The name becomes the ISSM marshalling filename and execution directory, so
    % it MUST be unique per PROCESS, not per job. `randi` is not enough: like
    % calc_steady_state, callers seed rng from SLURM_JOB_ID for reproducibility,
    % so two continuations launched inside ONE allocation draw the SAME number,
    % land on the same <name>.bin, and clobber each other mid-marshal. Observed
    % 2026-07-28: two legs launched together in job 36268159 both produced
    % `20260728_2152relaxleg2_rndi_52958.bin` and both died with
    % "CheckFile error message: error while looking in binary file. Found a
    % string of size 0". Nothing was corrupted (the crash precedes any write and
    % the lock released cleanly) but both legs were lost.
    %
    % feature('getpid') is distinct for concurrent MATLABs on one node; the
    % randi is kept so that sequential runs from the SAME pid still differ.
    md.miscellaneous.name = sprintf('%srelaxleg%d_pid%d_rndi_%d', ...
                                    char(date_str), leg_number, ...
                                    feature('getpid'), randi(65000));
    md.cluster = generic('name', oshostname(), ...
                         'executionpath', [exec_root filesep], ...
                         'np', num_proc);

    fprintf(['[relax_continue] leg %d: t %.4f -> %.4f (%g yr), dt %g, ' ...
             'output_frequency %g, np %d\n'], leg_number, start_time, ...
            md.timestepping.final_time, add_years, md.timestepping.time_step, ...
            md.settings.output_frequency, num_proc);

    % ----------------------------------------------------------------- 11. solve --
    md = solve(md, 'transient');
    new = md.results.TransientSolution;
    n_new = numel(new);
    fprintf('[relax_continue] leg produced %d snapshots, t = %.4f -> %.4f\n', ...
            n_new, double(new(1).time), double(new(end).time));

    % ------------------------------------------- 12. post-solve assertions -----

    % G12 -- the mesh-verbatim invariant. If this ever fires, something in the
    % solve remeshed and every IceMelt run seeded off the result would be at a
    % different resolution.
    if md.mesh.numberofvertices ~= nvert_in || ~isequal(md.mesh.x, mesh_x_in)
        error(['ISSMSteadyStateUtils.relax_continue [G12]: the mesh CHANGED ' ...
               'during the leg (%d -> %d vertices). A continuation must be ' ...
               'mesh-verbatim; refusing to write.'], nvert_in, ...
              md.mesh.numberofvertices);
    end

    % G13 -- the gate normalises MS_diff_H by mean(smb). A change here would
    % silently rescale the WHOLE historical series, not just the new leg.
    smb_mean_out = mean(md.smb.mass_balance);
    if abs(smb_mean_out - smb_mean_in) > 1e-12 * max(1, abs(smb_mean_in))
        error(['ISSMSteadyStateUtils.relax_continue [G13]: mean(smb.mass_balance) ' ...
               'changed %.12g -> %.12g. steady_state_gate.py normalises by it, ' ...
               'so this would rescale every historical point too.'], ...
              smb_mean_in, smb_mean_out);
    end

    % G10 -- the seam. A duplicated time reads to the gate as dt = 0 (a division
    % by zero) and a non-monotonic one as a spurious MS_diff_H spike.
    t_hist_end = double(history(end).time);
    t_new_1    = double(new(1).time);
    if ~(t_new_1 > t_hist_end)
        error(['ISSMSteadyStateUtils.relax_continue [G10]: seam violation -- ' ...
               'the new leg''s first snapshot is at t = %.6f but the history ' ...
               'ends at t = %.6f. The gate would see dt <= 0 there.'], ...
              t_new_1, t_hist_end);
    end

    % ---------------------------------------------------- 13. merge and thin ----
    merged = ISSMMiscUtils.concat_transientsol(history, new);
    t_all = arrayfun(@(s) double(s.time), merged);
    if any(diff(t_all) <= 0)
        bad = find(diff(t_all) <= 0, 1);
        error(['ISSMSteadyStateUtils.relax_continue [G10]: merged times are not ' ...
               'strictly increasing at index %d (%.6f -> %.6f).'], ...
              bad, t_all(bad), t_all(bad+1));
    end

    [retained, thin_info] = ISSMMiscUtils.thin_transientsol(merged, policy);
    fprintf('[relax_continue] retention: %s\n', thin_info.reason);

    md.results = struct('TransientSolution', retained);

    % --------------------------------------------------------- 14. provenance ---
    leg = struct( ...
        'leg',                      leg_number, ...
        'kind',                     'relax_continue', ...
        'start_time',               start_time, ...
        'final_time',               double(retained(end).time), ...
        'years',                    add_years, ...
        'time_step',                md.timestepping.time_step, ...
        'output_frequency',         md.settings.output_frequency, ...
        'amr_frequency',            md.transient.amr_frequency, ...
        'requested_outputs_sha256', sha256_of(strjoin(cellstr(md.transient.requested_outputs), '|')), ...
        'n_snapshots_solved',       n_new, ...
        'n_snapshots_retained',     numel(retained), ...
        'n_snapshots_history_in',   n_hist, ...
        'retention_thinned',        thin_info.thinned, ...
        'retention_reason',         thin_info.reason, ...
        'slurm_job_id',             getenv('SLURM_JOB_ID'), ...
        'host',                     oshostname(), ...
        'wallclock_s',              toc(t_wall_start), ...
        'helpers_git_sha',          helpers_sha, ...
        'helpers_dirty',            helpers_dirty, ...
        'issm_version',             num2str(issmversion()), ...
        'matlab_version',           version(), ...
        'num_proc',                 num_proc, ...
        'smb_mean',                 smb_mean_out, ...
        'n_vertices',               md.mesh.numberofvertices, ...
        'config_path',              config_path, ...
        'reconstructed',            false);

    prov.legs{end+1}    = leg;
    prov.retention      = policy;
    prov.cumulative     = recompute_cumulative(prov.legs);
    md.miscellaneous.dummy.relax_provenance = prov;

    fprintf(['[relax_continue] cumulative: %d legs, %g total relaxed yr ' ...
             '(%g mesh-fixed), t_final %.4f\n'], ...
            prov.cumulative.n_legs, prov.cumulative.total_relaxed_years, ...
            prov.cumulative.mesh_fixed_relaxed_years, ...
            prov.cumulative.current_final_time);

    % ------------------------------------------------------------- 15. write ---
    ISSMMiscUtils.save_model_atomic(md, run_dir, 'final_md.mat', keep_generations);

    % --------------------------------------------------- 16. sidecar + unlock ---
    % The .mat copy is authoritative; this is the greppable mirror so the history
    % can be read from shell without loading a 20 MB MAT or starting MATLAB.
    sidecar = fullfile(run_dir, 'relax_provenance.json');
    fid = fopen(sidecar, 'w');
    if fid > 0
        fprintf(fid, '%s\n', jsonencode(prov, 'PrettyPrint', true));
        fclose(fid);
        fprintf('[relax_continue] provenance sidecar: %s\n', sidecar);
    else
        fprintf(2, '[relax_continue] WARNING: could not write %s\n', sidecar);
    end

    provenance = prov;
    clear cleanup   % releases the lock
    fprintf('[relax_continue] DONE in %.1f s\n', toc(t_wall_start));
end


% =========================================================================== %
% locking (G4)
% =========================================================================== %

function acquire_lock(lock_path, config)
% Two concurrent legs racing on one in-place file would interleave a rename with
% a save and lose the run. The lock is advisory (a file, not flock) because that
% is what survives a node reboot legibly: the operator can read it.
%
% A STALE lock -- one whose job id is no longer in `squeue --me` -- is refused
% rather than silently stolen, because the alternative is that a hung-but-alive
% job gets its file swapped underneath it.
    if exist(lock_path, 'file')
        txt = fileread(lock_path);
        tok = regexp(txt, 'job_id:\s*(\S+)', 'tokens', 'once');
        held_by = '';
        if ~isempty(tok); held_by = tok{1}; end
        alive = job_is_alive(held_by);
        forced = isfield(config, 'FORCE') && config.FORCE;
        if alive
            error(['ISSMSteadyStateUtils.relax_continue [G4]: %s is held by ' ...
                   'job %s, which is STILL RUNNING. Two legs must not race on ' ...
                   'one in-place file. --force does NOT override a live lock.' ...
                   '\n--- lock ---\n%s'], lock_path, held_by, txt);
        end
        if ~forced
            error(['ISSMSteadyStateUtils.relax_continue [G4]: %s is present but ' ...
                   'its holder (job %s) is not in the queue -- a STALE lock ' ...
                   'from a crashed leg. Check that final_md.mat is intact ' ...
                   '(relax.py --check <RUN>) and re-run with --force.' ...
                   '\n--- lock ---\n%s'], lock_path, held_by, txt);
        end
        fprintf(2, ['[relax_continue] --force: breaking stale lock held by ' ...
                    'job %s\n'], held_by);
        delete(lock_path);
    end
    fid = fopen(lock_path, 'w');
    if fid < 0
        error('ISSMSteadyStateUtils.relax_continue [G4]: cannot create %s', lock_path);
    end
    fprintf(fid, 'job_id: %s\nhost: %s\npid: %d\nstarted: %s\n', ...
            getenv('SLURM_JOB_ID'), oshostname(), feature('getpid'), ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fclose(fid);
end

function release_lock(lock_path)
    if exist(lock_path, 'file')
        delete(lock_path);
    end
end

function alive = job_is_alive(job_id)
% Empty job id means the previous holder ran outside Slurm; treat as dead, since
% there is no way to check and the stale-lock branch already requires --force.
    alive = false;
    if isempty(job_id) || strcmp(job_id, '<none>')
        return
    end
    [rc, out] = system(sprintf('squeue -h -j %s -o %%i 2>/dev/null', job_id));
    alive = (rc == 0) && ~isempty(strtrim(out));
end


% =========================================================================== %
% provenance helpers
% =========================================================================== %

function prov = existing_provenance(md, run_dir)
% Read the provenance already inside the .mat, or reconstruct a leg-1 record if
% this is the first continuation of a run that predates provenance.
%
% THE TRAP THIS MAKES VISIBLE (design §8): the old flatbed seed has
% start_time = 2625 / final_time = 2700, i.e. it was ITSELF a hand-made
% continuation chunk. The reconstructed record therefore reports 75 yr, which is
% honest -- it is what the file can PROVE. The ~1700 yr figure in the handoff is
% a claim ABOUT the file, not a claim FROM it. Nothing here guesses.
    if isstruct(md.miscellaneous.dummy) && ...
            isfield(md.miscellaneous.dummy, 'relax_provenance')
        prov = md.miscellaneous.dummy.relax_provenance;
        if ~iscell(prov.legs)
            % jsondecode / a round-trip can turn a 1x1 cell into a bare struct.
            prov.legs = num2cell(prov.legs);
        end
        fprintf('[relax_continue] found existing provenance: %d leg(s)\n', ...
                numel(prov.legs));
        return
    end

    sols = md.results.TransientSolution;
    leg1 = struct( ...
        'leg', 1, 'kind', 'calc_steady_state:final', ...
        'start_time', double(sols(1).time), ...
        'final_time', double(sols(end).time), ...
        'years', double(sols(end).time) - double(sols(1).time), ...
        'time_step', md.timestepping.time_step, ...
        'output_frequency', md.settings.output_frequency, ...
        'amr_frequency', md.transient.amr_frequency, ...
        'requested_outputs_sha256', [], ...
        'n_snapshots_solved', numel(sols), ...
        'n_snapshots_retained', numel(sols), ...
        'n_snapshots_history_in', 0, ...
        'retention_thinned', false, ...
        'retention_reason', 'predates the retention policy', ...
        'slurm_job_id', [], 'host', [], 'wallclock_s', [], ...
        'helpers_git_sha', [], 'helpers_dirty', [], ...
        'issm_version', [], 'matlab_version', [], 'num_proc', [], ...
        'smb_mean', mean(md.smb.mass_balance), ...
        'n_vertices', md.mesh.numberofvertices, ...
        'config_path', [], ...
        'reconstructed', true);
    prov = struct('schema_version', 1, ...
                  'run_dir', run_dir, ...
                  'origin', struct('generator', 'calc_steady_state', ...
                                   'config_path', [], 'config_sha256', []), ...
                  'cumulative', struct(), ...
                  'retention', [], ...
                  'legs', {{leg1}});
    prov.cumulative = recompute_cumulative(prov.legs);
    fprintf(['[relax_continue] no provenance in the .mat; RECONSTRUCTED leg 1 ' ...
             'from the file (t %.4f -> %.4f, %g yr, %d snapshots). Marked ' ...
             'reconstructed = true. Anything not knowable from the file (the ' ...
             'Slurm job, the helpers SHA) is left null, never guessed.\n'], ...
            leg1.start_time, leg1.final_time, leg1.years, leg1.n_snapshots_solved);
end

function c = recompute_cumulative(legs)
% Recomputed from the leg list every time, never incremented, so it cannot drift
% away from the legs that justify it.
%
% mesh_fixed_relaxed_years is kept SEPARATE from total_relaxed_years because the
% AMR transient and the fix-icefront relax are not the same kind of year as a
% mesh-fixed relax -- conflating them is precisely how "250 yr" became ambiguous.
    years = 0; mesh_fixed = 0; starts = []; finals = []; recon = false;
    for k = 1:numel(legs)
        L = legs{k};
        years = years + double(L.years);
        if isfield(L, 'amr_frequency') && ~isempty(L.amr_frequency) && L.amr_frequency == 0
            mesh_fixed = mesh_fixed + double(L.years);
        end
        starts(end+1) = double(L.start_time); %#ok<AGROW>
        finals(end+1) = double(L.final_time); %#ok<AGROW>
        recon = recon || (isfield(L, 'reconstructed') && ~isempty(L.reconstructed) && L.reconstructed);
    end
    c = struct('n_legs', numel(legs), ...
               'total_relaxed_years', years, ...
               'mesh_fixed_relaxed_years', mesh_fixed, ...
               'first_start_time', min(starts), ...
               'current_final_time', max(finals), ...
               'reconstructed', recon);
end

function [sha, dirty] = helpers_git_state()
    here = fileparts(fileparts(mfilename('fullpath')));   % .../ISSMHelperFunctions
    sha = ''; dirty = false;
    [rc, out] = system(sprintf('cd %s && git rev-parse --short HEAD 2>/dev/null', here));
    if rc == 0; sha = strtrim(out); end
    [rc, out] = system(sprintf('cd %s && git status --porcelain 2>/dev/null', here));
    if rc == 0; dirty = ~isempty(strtrim(out)); end
end

function h = sha256_of(str)
% Hash of the requested_outputs list, so a leg that silently ran a DIFFERENT
% output list is detectable after the fact even though nothing here changes it.
    try
        md5 = java.security.MessageDigest.getInstance('SHA-256');
        b = md5.digest(uint8(str));
        h = lower(reshape(dec2hex(typecast(b, 'uint8'))', 1, []));
    catch
        h = [];
    end
end
