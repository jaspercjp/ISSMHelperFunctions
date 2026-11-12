%CALC_STEADY_STATE simulates a steady state ice sheet / ice shelf given a configuration
% file found at config_dir. The configuration file contains all necessary information needed
% to set up the simulation and run it.
function [md, config] = calc_steady_state(config_dir, varargin)
    % Use this script to calculate a steady state with a fixed-slope bed.
    % Note that here we are using the 'RestartModel' argument of RunTransient even
    % though we are not really restarting a model.
    p = inputParser;
    addParameter(p, 'CoarseSteadyModel', [], @(x) isempty(x) || isa(x, 'model'));
    addParameter(p, 'RunAMRInterp', true    , @(x) islogical(x));
    addParameter(p, 'IsThermal', false, @(x) islogical(x));
    parse(p,varargin{:});
    
    COARSE_STEADY_MODEL  = p.Results.CoarseSteadyModel;
    RUN_AMR_INTERP       = p.Results.RunAMRInterp;
    IS_THERMAL           = p.Results.IsThermal;

    % Ensure true pseudorandomness in file ID when running on slurm
    jobID_str = getenv('SLURM_JOB_ID');
    if isempty(jobID_str)
        rng('shuffle');
    else
        jobID_num = str2double(jobID_str);
        rng(jobID_num);
    end

    config=jsondecode(fileread(config_dir));
    % Output + execution roots come from the environment (exported by the
    % launcher: $DATA_ROOT / $ISSM_EXECUTION_ROOT), with a legacy cluster
    % fallback. Steady states are saved under <DATA_ROOT>/steady_states/.
    data_root = getenv('DATA_ROOT');
    if isempty(data_root)
        error('ISSMHelperFunctions:noDataRoot', 'DATA_ROOT is not set; the launchers (Core.issm_run, _calc_steady_state.sh) export it.');
    end
    exec_root = getenv('ISSM_EXECUTION_ROOT');
    if isempty(exec_root)
        error('ISSMHelperFunctions:noExecRoot', 'ISSM_EXECUTION_ROOT is not set; the launchers (Core.issm_run, _calc_steady_state.sh) export it.');
    end
    config.SAVE_PATH = [fullfile(data_root, 'steady_states') filesep];
    config.EXECUTION_PATH = [exec_root filesep];
    config.SAVE_PATH=char(config.SAVE_PATH);
    config.EXECUTION_PATH=char(config.EXECUTION_PATH);
    % Domain outline (.exp): referenced relative to the config file's folder
    % (e.g. "domains/Outline_...exp"), or as an absolute path.
    df = char(config.DOMAIN_FILE);
    if startsWith(df, '/')
        TOTAL_DOMAIN_FILEPATH = df;
    else
        TOTAL_DOMAIN_FILEPATH = fullfile(fileparts(config_dir), df);
    end
    % Hand the already-resolved absolute outline path downstream: interp_refine_md3d
    % re-meshes from the same outline and has no way to know config_dir otherwise.
    config.TOTAL_DOMAIN_FILEPATH = char(TOTAL_DOMAIN_FILEPATH);

    date_str = ISSMMiscUtils.get_date_str();

    domain_str_split = split(config.DOMAIN_FILE, '/');
    domain_str = domain_str_split{end};
    domain_str_split = split(strip(domain_str),'.');
    domain_str = domain_str_split{1};
    domain_str_split = split(domain_str, '_');
    width_str = domain_str_split{end};

    % Two DISTINCT mesh-size specs live in this config, for two different
    % purposes -- do not let one silently clobber the other (see the fix for
    % the flatbed_L320km_W80km_smb60cmyr crash, "amr.icefront_resolution
    % should have values below 5000"):
    %   - config.HMAX/HMIN        : bounds for the initial coarse bamg mesh
    %                               generation, used once, right here.
    %   - config.AMR_SETTING.hmax/hmin : bounds for the AMR *refinement*
    %                               phase (md.amr.hmax/hmin), read once into
    %                               AMR_HMAX/AMR_HMIN below and used wherever
    %                               md.amr.hmax/hmin needs to be (re)assigned.
    % Each is read from the config exactly once, into a locally-scoped
    % variable, and never overwritten from the other source.
    COARSE_MESH_HMAX = config.HMAX;
    COARSE_MESH_HMIN = config.HMIN;

    if ~isa(COARSE_STEADY_MODEL, 'model')
        md=bamg(model, 'domain', TOTAL_DOMAIN_FILEPATH, 'hmax', COARSE_MESH_HMAX, 'hmin', COARSE_MESH_HMIN);
    else
        md=COARSE_STEADY_MODEL;
    end

    if strcmp(config.PARAM_TYPE, 'ConstSlope')
        PARAM_FUNC = @(m) ISSMParamUtils.ParamConstSlope(m, config.PARAM_CONFIG.SLOPE, ...
                                                            'Config', config, ...
                                                            'InitGeometry', true, ...
                                                            'Offset', config.PARAM_CONFIG.OFFSET);
        if IS_THERMAL
            FOLDER_NAME = [char(date_str) sprintf('_CONSTSLOPE_steady_%dyr_smb%dcmyr_slope%de-4_%s_surfTn%d/', ... 
                            config.NUM_YRS,... 
                            config.SURFACE_MASS_BALANCE*100, ...
                            config.PARAM_CONFIG.SLOPE*1e4, ...
                            width_str, ...
                            abs(config.THERMAL_3D_SETTINGS.surface_temp))];
        else
            FOLDER_NAME = [char(date_str) sprintf('_CONSTSLOPE_steady_%dyr_smb%dcmyr_slope%de-4_%s/', ... 
                            config.NUM_YRS,... 
                            config.SURFACE_MASS_BALANCE*100, ...
                            config.PARAM_CONFIG.SLOPE*1e4, ...
                            width_str)];
        end
    elseif strcmp(config.PARAM_TYPE, 'MISMIP')
        if isfield(config, 'STARTUP_MODEL')
            startup_load_path = config.SAVE_PATH;   % prior steady states live here
            startup_data = load([startup_load_path config.STARTUP_MODEL]);
            startup_model = startup_data.md;
        else
            startup_model = [];
        end
        % `config` (not config.PARAM_CONFIG) is passed through so that top-level
        % keys such as FRICTION_LAW reach ParamMISMIP.
        PARAM_FUNC = @(m) ISSMParamUtils.ParamMISMIP(m, 'Config', config, 'NoSlip', true, 'StartupModel', startup_model);
        if isfield(config, 'FRICTION_LAW')
            friction_law = config.FRICTION_LAW;
        else
            disp('[calc_steady_state] :: No friction parametrization specified. Defaulting to Weertman.')
            friction_law = 'frictionweertman';
        end
        FOLDER_NAME = [char(date_str) sprintf('_MISMIP_%s_steady_smb%dcmyr/',...
                                               friction_law, ...
                                               config.SURFACE_MASS_BALANCE*100)];
    else
        error(['The parametrization ' config.PARAM_TYPE ' is not currently supported!'])
    end
    md=PARAM_FUNC(md);
    
    % If provided, initialize geometry from the 1D steady state solution
    if isfield(config, 'INIT_STEADY_1D_FILE')
        % 1D width-averaged seed solutions live in the repo under
        % SteadyState/1D_steady_sol/, i.e. ../../1D_steady_sol relative to the
        % STEADY_STATE_CONFIGS/ folder holding this config file.
        total_1D_filepath = fullfile(fileparts(config_dir), '..', '..', '1D_steady_sol', config.INIT_STEADY_1D_FILE);
        file_exists = exist(total_1D_filepath, 'file');
        if file_exists
            % SEED_Y_AWARE (default false = historical y-uniform behaviour).
            % When true the 1-D profile is mapped onto the 2-D bed by preserving
            % its SURFACE over grounded ice instead of its thickness, which
            % removes the spurious transverse surface step a laterally varying
            % bed otherwise imprints (+499.5 m on MISMIP+; see
            % init_steady_from_1Dwavg.m's header). Only meaningful where the bed
            % varies in y -- on flatbed/constslope the two paths coincide.
            if isfield(config, 'SEED_Y_AWARE') && config.SEED_Y_AWARE
                seed_bed = md.geometry.bed;
            else
                seed_bed = [];
            end
            [h, s, b, ~] = ISSMSteadyStateUtils.init_steady_from_1Dwavg(md.mesh.x, md.mesh.y, total_1D_filepath, 1, seed_bed);
            md.geometry.thickness = h;
            md.geometry.base = b;
            md.geometry.surface = s;
            md = ISSMModelUtils.adjust_thk_ols(md);
        else
            warning(['The specified 1D steady state solution ' total_1D_filepath ' was not found']);
        end
    end

    % Generate a coarse steady state model 
    if ~isa(COARSE_STEADY_MODEL, 'model')
        md=setflowequation(md,'SSA','all');
        md.masstransport.stabilization = 1;
        md.levelset.stabilization      = 1;
        md.timestepping.start_time = 0;
        md.timestepping.final_time = config.NUM_YRS;
        md.timestepping.time_step  = config.TIME_STEP; 
        md.settings.output_frequency = config.OUTPUT_FREQ;
        % Optional per-config override of ISSM's default solver_residue_threshold
        % (1e-6). Added 2026-07-28 for MISMIP_L640km_W80km.json (Weertman): the
        % coarse transient hits a real but self-limiting stiffness a few
        % timesteps in (residue ~1e-5-1e-4 around t=2-4yr, both with and
        % without INIT_STEADY_1D_FILE), which the default 1e-6 bar rejects even
        % though the solve recovers and runs cleanly for the following ~95 yr
        % once past it (verified: 500-iteration/100yr coarse-only run at 1e-4
        % with zero further residue violations after iteration ~22). Only
        % covers THIS (coarse) solve -- the AMR-refinement/relax solves further
        % down are not covered by this override and would need it added
        % separately if the same instability turns up there.
        if isfield(config, 'SOLVER_RESIDUE_THRESHOLD')
            md.settings.solver_residue_threshold = config.SOLVER_RESIDUE_THRESHOLD;
        end
        md.miscellaneous.name = [char(date_str) 'coarse_1D_rndi_' num2str(randi(65000))];
        md.cluster = generic('name', oshostname, 'executionpath', config.EXECUTION_PATH, 'np', config.NUM_PROC);
        md = solve(md,'transient');
        TOTAL_SAVE_PATH = [config.SAVE_PATH FOLDER_NAME '/'];
        if ~exist(TOTAL_SAVE_PATH, 'dir')
            system(['mkdir -p ' TOTAL_SAVE_PATH]);
        end
        save([TOTAL_SAVE_PATH 'coarse_md.mat'], 'md');
    end

    if IS_THERMAL
        md_3d = ISSMModelUtils.ExtrudeRunThermal(md, config, 'FolderName', FOLDER_NAME);
        if RUN_AMR_INTERP
            md_3d = ISSMSteadyStateUtils.interp_refine_md3d(md_3d, config, 'FolderName', FOLDER_NAME);
        end
    elseif RUN_AMR_INTERP
        % The case of 2D steady state (no heat transport), interpolate with AMR 
        disp('[calc_steady_state] :: Coarse steady state genereated, proceeding to use AMR to refine...')
        md.miscellaneous.name = [char(date_str) 'steady_amrinterp_rdni_' num2str(randi(65000))];
        md = transientrestart(md);
        md.timestepping.final_time = md.timestepping.final_time + config.AMR_INTERP_SETTINGS.NUM_YRS;
        md.transient.amr_frequency = config.AMR_SETTING.amr_freq;
        md.settings.output_frequency    = config.AMR_INTERP_SETTINGS.OUTPUT_FREQ;
        % AMR refinement mesh-size bounds: read from config.AMR_SETTING exactly
        % once, into AMR_HMAX/AMR_HMIN, here. These bound md.amr.hmax/hmin for
        % BOTH AMR-driven solves below (this transient and the mesh-fixing
        % relax further down) -- never from config.HMAX/HMIN, which are a
        % separate spec for the one-time initial coarse mesh generation above.
        AMR_HMAX = config.AMR_SETTING.hmax;
        AMR_HMIN = config.AMR_SETTING.hmin;
        setting_name=fieldnames(config.AMR_SETTING);
        for i=1:length(setting_name)
            if (~strcmp(setting_name{i},'amr_freq')) && ~strcmp(setting_name{i},'hmin') && ~strcmp(setting_name{i},'hmax')
                disp(setting_name{i})
                md.amr.(setting_name{i}) = config.AMR_SETTING.(setting_name{i});
            end
        end
        % Assigned here, immediately before this phase's solve, from the
        % single AMR_SETTING-sourced variables above.
        md.amr.hmax = AMR_HMAX;
        md.amr.hmin = AMR_HMIN;
        md.amr.restart = 1;
        req_outputs = {'default', 'StrainRatexx', 'StrainRateyy', ...
                   'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy', ...
                   'MaterialsRheologyB', 'Bed', 'FrictionAlpha2'};
        md.transient.requested_outputs=req_outputs;
        md = solve(md, 'transient');
        md = ISSMMiscUtils.delete_prev_transientsol(md);
        % SAVE_AMR_MD (added 2026-07-28, Docs/Design-relax-continuation.md §10.2).
        % OPTIONAL; ABSENT == 'full' == exactly the behaviour that shipped before
        % this key existed, so a job that re-reads this file mid-run is
        % unaffected.
        %
        %   'full' : today's behaviour -- every AMR snapshot. For debugging an
        %            AMR problem. 400-800 MB per run.
        %   'thin' : the model plus the two snapshots amr_md.mat is ever actually
        %            used for (the last one, and the last one carrying a mesh --
        %            see ISSMMiscUtils.thin_amr_snapshots). ~10-20 MB, with
        %            resumability preserved in full.
        %   'none' : skip the write entirely.
        %
        % Nothing on the live path reads this file from disk: the AMR mesh
        % reaches the next phase IN MEMORY via copy_amr_mesh (see
        % `amrinterp_endsol` two lines below, which is taken from `md`, NOT from
        % the file), the gate refuses to grade the AMR phase, and
        % compare_steady_states loads only coarse and final. The thinning is
        % applied to a COPY so `md` -- and therefore copy_amr_mesh -- is untouched.
        AMR_SAVE_MODE = 'full';
        if isfield(config, 'SAVE_AMR_MD') && ~isempty(config.SAVE_AMR_MD)
            AMR_SAVE_MODE = char(config.SAVE_AMR_MD);
        end
        switch lower(AMR_SAVE_MODE)
            case 'none'
                disp('[calc_steady_state] :: SAVE_AMR_MD=none, skipping amr_md.mat');
            case 'thin'
                % The helper takes `md` BY VALUE, so MATLAB's copy-on-write means
                % thinning inside it cannot touch the `md` in this scope -- which
                % matters, because `amrinterp_endsol` below needs the full array.
                % It also saves under the variable name `md`, so every existing
                % reader's `load(..., 'md')` keeps working.
                ISSMMiscUtils.save_thin_amr_md(TOTAL_SAVE_PATH, md);
            otherwise
                save([TOTAL_SAVE_PATH 'amr_md.mat'], 'md');
        end
        % Whole array, not (end): only snapshots coinciding with an AMR remesh
        % event carry MeshX/MeshY/MeshElements, so copy_amr_mesh scans backwards
        % for the last one that does. Passing (end) crashed job 36121382.
        amrinterp_endsol = md.results.TransientSolution;

        % Relax the model by uniformly refining mesh one last time
        md.miscellaneous.name = [char(date_str) 'fixicefrontrelax_rdni_' num2str(randi(65000))];
        md=ISSMModelUtils.ReinitializeAMRModel(md, 'ParamFunc', PARAM_FUNC);
        md=ISSMParamUtils.init_BC(md, 'NoSlip', true);
        md.timestepping.start_time = md.timestepping.final_time; 
        % AMR_INTERP_SETTINGS.NUM_YRS is the TOTAL post-coarse relax, split into
        % the two /2 legs here and below (500 -> 2 x 250 yr for the standardised
        % configs). The old comment said "100 years", which was only ever true
        % for a NUM_YRS of 200.
        md.timestepping.final_time = md.timestepping.final_time + config.AMR_INTERP_SETTINGS.NUM_YRS/2;
        md = ISSMMeshUtils.copy_amr_mesh(md, amrinterp_endsol);
        md.transient.amr_frequency = 100000;
        % Re-assigned here, immediately before this second AMR-driven solve,
        % from the SAME AMR_HMAX/AMR_HMIN read once above -- not re-read from
        % config, and not from config.HMAX/HMIN.
        md.amr.hmin = AMR_HMIN; md.amr.hmax = AMR_HMAX;
        md = solve(md,'transient');
        fixicefront_end_sol = md.results.TransientSolution;  % whole array; see above

        % Relax model fixing the mesh
        md=ISSMModelUtils.ReinitializeAMRModel(md, 'ParamFunc', PARAM_FUNC);
        md=ISSMParamUtils.init_BC(md, 'NoSlip', true);
        md.timestepping.start_time = md.timestepping.final_time; 
        % Second (and last) /2 leg -- see the note above. This is the leg that
        % becomes final_md: delete_prev_transientsol below drops the earlier
        % legs, which is why steady_state_gate.py's "final" phase spans only
        % NUM_YRS/2 years (5 -> 250) and not the full NUM_YRS.
        md.timestepping.final_time = md.timestepping.final_time + config.AMR_INTERP_SETTINGS.NUM_YRS/2;
        md = ISSMMeshUtils.copy_amr_mesh(md, fixicefront_end_sol);
        md.transient.amr_frequency = 0;
        md = solve(md,'transient');
        md = ISSMMiscUtils.delete_prev_transientsol(md);
        % RETENTION_FILE (added 2026-07-28, Docs/Design-relax-continuation.md
        % §4.6). OPTIONAL; ABSENT == no thinning == exactly the behaviour that
        % shipped before this key existed, so a job that re-reads this file
        % mid-run is unaffected.
        %
        % WHY FRESH RUNS SHOULD USE THE SAME POLICY AS CONTINUATIONS. Without
        % this, leg 1 is 192 MB and leg 2 is 11 MB -- an absurd discontinuity,
        % and it means the FIRST continuation silently destroys ~180 MB of
        % history that the policy would have kept. One kernel, two call sites.
        %
        % The policy is a PATH, not inline JSON, so there is exactly one
        % materialisation of it in the repo
        % (Scripts/Workflows/SteadyState/Configs/retention_standard.json, itself
        % pinned against Core.relax_retention.STANDARD_POLICY by
        % Scripts/Core/Tests/test_relax_retention.py). Relative paths resolve
        % against the config file's own folder, like DOMAIN_FILE.
        if isfield(config, 'RETENTION_FILE') && ~isempty(config.RETENTION_FILE)
            rf = char(config.RETENTION_FILE);
            if ~startsWith(rf, '/')
                rf = fullfile(fileparts(config_dir), rf);
            end
            if ~exist(rf, 'file')
                error(['[calc_steady_state] RETENTION_FILE %s not found. ' ...
                       'Refusing to fall back on "no thinning" silently -- ' ...
                       'that would produce a 192 MB leg-1 file that the first ' ...
                       'continuation then thins away.'], rf);
            end
            retention_policy = jsondecode(fileread(rf));
            [thinned_sols, thin_info] = ISSMMiscUtils.thin_transientsol( ...
                md.results.TransientSolution, retention_policy);
            disp(['[calc_steady_state] :: retention: ' thin_info.reason]);
            md.results.TransientSolution = thinned_sols;
        end
        save([TOTAL_SAVE_PATH 'final_md.mat'], 'md');
    end
    
    disp('[calc_steady_state] :: Steady state generated!')
end