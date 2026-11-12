function [out, info] = thin_transientsol(sols, policy)
%THIN_TRANSIENTSOL  Apply the relax-continuation retention policy to a
%   TransientSolution struct array.
%
%   [OUT, INFO] = ISSMMiscUtils.thin_transientsol(SOLS, POLICY)
%
%   SOLS   1xN struct array, ordered oldest -> newest (ISSM's own ordering).
%   POLICY struct decoded from the run config's RETENTION block. It MUST carry
%          seed_snapshots, keep_fields, window_yr, window_cadence_yr,
%          archive_cadence_yr and min_window_points. There are DELIBERATELY no
%          defaults here -- see "WHY NO DEFAULTS" below.
%
%   OUT    the retained struct array, same field NAMES as SOLS.
%   INFO   struct: .thinned .reason .n_in .n_out .n_seed .n_thin
%                  .window_points .kept_idx
%
% ---------------------------------------------------------------------------
% WHY NO DEFAULTS
% ---------------------------------------------------------------------------
% The policy object's single source of truth is Core.relax_retention.STANDARD_POLICY
% (Python). relax.py bakes it verbatim into the derived config JSON that this
% function's caller reads. If this file carried its own copy of the numbers, the
% two would drift -- this repo already has that scar (Scripts/Core/issm_env.sh
% exists because one shell recipe lived in two files and diverged inside a day).
% So a missing field is a hard ERROR, not a fallback.
%
% ---------------------------------------------------------------------------
% THE SELECTION RULE (mirrored, and unit-tested, in
% Core.relax_retention.plan_retention / Scripts/Core/Tests/test_relax_retention.py)
% ---------------------------------------------------------------------------
%   1. The last `seed_snapshots` records are kept COMPLETE. transientrestart --
%      and therefore every IceMelt experiment that seeds off this run -- reads
%      md.results.TransientSolution(end) and nothing else, so thinning the seed
%      would install empty velocity/geometry fields and the next leg would start
%      from garbage. Guard G2 in relax_continue.m asserts this did not happen.
%   2. Walk BACKWARDS from the seed. Keep a candidate when its time is at least
%      one cadence older than the last record already kept, where the cadence is
%      window_cadence_yr inside the trailing window_yr and archive_cadence_yr
%      outside it. Backwards, not forwards, so the cadence is anchored at the
%      NEWEST end -- the end the gate's 300 yr trend fit actually uses. A forward
%      walk would dump the ragged remainder there instead.
%   3. Retained non-seed records keep only `keep_fields`; every other field is
%      set to [] (see "THE MATLAB CONSTRAINT" below).
%   4. G6: if fewer than min_window_points retained records land inside the
%      trailing window_yr, DISCARD THE PLAN and return SOLS unchanged. A file
%      steady_state_gate.py cannot grade is worse than a file that is too big:
%      the gate needs > 2 points to fit a trend at all, and with no trend it
%      cannot detect a RISING series, so it would pass a diverging run on level
%      alone. The operator can always thin later; they cannot un-thin.
%
% ---------------------------------------------------------------------------
% THE MATLAB CONSTRAINT, AND THE FOOTGUN IT LEAVES BEHIND
% ---------------------------------------------------------------------------
% MATLAB struct arrays require IDENTICAL field sets across elements, so a "full"
% record and a "thin" record cannot be concatenated as-is. Two ways out were
% considered (design §4.5):
%
%   (a) uniform field set with [] in the dropped fields  <-- CHOSEN
%   (b) a separate md.results.RelaxHistory array
%
% (a) wins because it leaves FOUR consumers working unmodified: steady_state_gate.py,
% compare_steady_states.py, check_steady_state.m and ISSM's own transientrestart.
% (b) is semantically cleaner but would need all four updated and would give
% every archived run a different shape from every new one.
%
% The cost of (a), stated plainly: isfield(sol,'Vx') is TRUE on a thin record even
% though Vx is empty. So transientrestart(md, step) applied to a THIN step would
% silently install md.initialization.vx = []. Guard G2 makes the tail record
% provably full, and nothing in this repo passes an explicit `step`, but the
% footgun is real and this comment is the warning.
%
% VERIFIED 2026-07-28: a saved model carrying 0x0 doubles inside
% TransientSolution elements round-trips through Core.ModelIO.issm_mat --
% transient_solutions() returns them as empty arrays and does not raise, and
% steady_state_gate.ms_diff_h computes cleanly over the result. (Design §12
% listed this as unverified.)

    % ---------------------------------------------------------------- guards --
    if isempty(sols)
        out = sols;
        info = struct('thinned', false, 'reason', 'empty TransientSolution', ...
                      'n_in', 0, 'n_out', 0, 'n_seed', 0, 'n_thin', 0, ...
                      'window_points', 0, 'kept_idx', []);
        return
    end
    required = {'seed_snapshots', 'window_yr', 'window_cadence_yr', ...
                'archive_cadence_yr', 'min_window_points'};
    for k = 1:numel(required)
        if ~isfield(policy, required{k})
            error(['ISSMMiscUtils.thin_transientsol: the RETENTION policy is ' ...
                   'missing required field "%s". This function carries NO ' ...
                   'defaults on purpose -- the policy''s single source of truth ' ...
                   'is Core.relax_retention.STANDARD_POLICY, which relax.py ' ...
                   'bakes into the derived config. A partial policy means the ' ...
                   'config was hand-written; fix the config.'], required{k});
        end
    end
    if policy.seed_snapshots < 1
        error(['ISSMMiscUtils.thin_transientsol: seed_snapshots = %d < 1. ' ...
               'transientrestart reads TransientSolution(end) and nothing ' ...
               'else, so a zero-seed policy destroys the run.'], ...
              policy.seed_snapshots);
    end

    n = numel(sols);
    t = zeros(1, n);
    for k = 1:n
        t(k) = double(sols(k).time);
    end
    if any(diff(t) <= 0)
        bad = find(diff(t) <= 0, 1);
        error(['ISSMMiscUtils.thin_transientsol: snapshot times are not ' ...
               'strictly increasing (t(%d) = %.6g, t(%d) = %.6g). Thinning a ' ...
               'non-monotonic series would produce a dt <= 0 in the gate''s ' ...
               'first difference, i.e. a divide-by-zero or a spurious spike.'], ...
              bad, t(bad), bad+1, t(bad+1));
    end

    % ------------------------------------------------------------ selection --
    n_seed      = min(double(policy.seed_snapshots), n);
    seed_idx    = (n - n_seed + 1):n;
    t_end       = t(end);
    window_start = t_end - double(policy.window_yr);
    win_cad     = double(policy.window_cadence_yr);
    arc_cad     = double(policy.archive_cadence_yr);

    thin_mask = false(1, n);
    last_kept_t = t(seed_idx(1));
    for i = (seed_idx(1) - 1):-1:1
        if t(i) > window_start
            cadence = win_cad;
        else
            cadence = arc_cad;
        end
        if cadence <= 0
            % cadence 0 means "no decimation in this tier" (the `full` profile).
            thin_mask(i) = true;
            last_kept_t  = t(i);
            continue
        end
        % 1e-9 relative slack: snapshot times are floating-point accumulations
        % of dt*output_frequency and land a hair under an exact multiple.
        if (last_kept_t - t(i)) >= cadence * (1 - 1e-9)
            thin_mask(i) = true;
            last_kept_t  = t(i);
        end
    end

    keep_mask            = thin_mask;
    keep_mask(seed_idx)  = true;
    kept_idx             = find(keep_mask);
    window_points        = sum(t(kept_idx) >= window_start);

    % ------------------------------------------------------------------ G6 ---
    if window_points < double(policy.min_window_points)
        out  = sols;
        info = struct( ...
            'thinned', false, ...
            'reason', sprintf(['G6 REFUSED to thin: only %d retained record(s) ' ...
                'would land in the trailing %g yr window, below ' ...
                'min_window_points = %d. steady_state_gate.py needs > 2 points ' ...
                'to fit a trend at all and cannot detect a RISING series ' ...
                'without one, so an under-populated window would silently pass ' ...
                'a diverging run. Keeping all %d snapshots unthinned instead.'], ...
                window_points, policy.window_yr, policy.min_window_points, n), ...
            'n_in', n, 'n_out', n, 'n_seed', n, 'n_thin', 0, ...
            'window_points', sum(t >= window_start), 'kept_idx', 1:n);
        warning('ISSMMiscUtils:thin_transientsol:G6', '%s', info.reason);
        return
    end

    % --------------------------------------------------------- field blanking --
    out = sols(kept_idx);
    if isfield(policy, 'keep_fields') && ~isempty(policy.keep_fields)
        keep_fields = policy.keep_fields;
        if ~iscell(keep_fields)
            keep_fields = cellstr(keep_fields);
        end
        % Structural bookkeeping fields are free (scalars/short strings) and
        % several consumers index by them, so they are never blanked.
        keep_fields = [keep_fields(:)', {'time', 'step', 'SolutionType'}];
        all_fields  = fieldnames(out);
        drop        = all_fields(~ismember(all_fields, keep_fields));
        % Blank rather than rmfield: a uniform field set is what makes a thin
        % record concatenable with a full one (see "THE MATLAB CONSTRAINT").
        is_seed = ismember(kept_idx, seed_idx);
        for k = 1:numel(out)
            if is_seed(k)
                continue    % the seed stays COMPLETE, bit for bit
            end
            for f = 1:numel(drop)
                out(k).(drop{f}) = [];
            end
        end
    end

    info = struct( ...
        'thinned', true, ...
        'reason', sprintf(['retained %d of %d snapshots (%d full seed + %d ' ...
            'thin); %d inside the trailing %g yr window (min %d)'], ...
            numel(kept_idx), n, n_seed, numel(kept_idx) - n_seed, ...
            window_points, policy.window_yr, policy.min_window_points), ...
        'n_in', n, 'n_out', numel(kept_idx), 'n_seed', n_seed, ...
        'n_thin', numel(kept_idx) - n_seed, ...
        'window_points', window_points, 'kept_idx', kept_idx);
end
