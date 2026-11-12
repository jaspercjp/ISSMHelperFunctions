function info = save_model_atomic(md, run_dir, filename, keep_generations)
%SAVE_MODEL_ATOMIC  Replace a model file in place without a window in which the
%   good copy is gone and the new one is not yet whole.
%
%   INFO = ISSMMiscUtils.save_model_atomic(MD, RUN_DIR, FILENAME, KEEP_GENERATIONS)
%
%   MD                the model to write
%   RUN_DIR           directory holding FILENAME
%   FILENAME          e.g. 'final_md.mat'
%   KEEP_GENERATIONS  how many previous copies to retain as
%                     final_md.leg<N>.mat (0 = none). Default 1.
%
%   INFO  struct: .path .tmp_path .generation_path .bytes .verified
%
% ---------------------------------------------------------------------------
% BE EXPLICIT ABOUT THE RISK
% ---------------------------------------------------------------------------
% final_md.mat is the product of 4-17 h of solve and is the seed every IceMelt
% experiment inherits VERBATIM. Overwriting it in place is destructive by
% construction. This function reduces, but does not eliminate, the exposure. It
% is the only thing standing between a crash and the loss of that compute, which
% is why it is this heavily commented.
%
% THE SEQUENCE
% ------------
%   1  tmp = RUN_DIR/<stem>.tmp-<SLURM_JOB_ID or pid>.mat
%      -- the trailing .mat is LOAD-BEARING, not cosmetic. MATLAB's load()
%      refuses any path whose extension is not .mat ("Input must be a MAT-file
%      or an ASCII file containing numeric data"), so a tmp named
%      final_md.mat.tmp-<id> writes fine and then fails the G14 verify reload.
%      Found the hard way, by G14, on 2026-07-28.
%   2  save(tmp, 'md')
%   3  reload tmp and verify it (G14)
%   4  if KEEP_GENERATIONS > 0: movefile(FILENAME -> <stem>.leg<N>.mat)
%   5  movefile(tmp -> FILENAME)
%   6  prune generations beyond KEEP_GENERATIONS, oldest first
%
% Steps 4 and 5 are renames WITHIN ONE DIRECTORY on ONE Lustre filesystem, so
% each is atomic: a concurrent reader sees either the whole old file or the whole
% new one, never a partial. That is the entire point. The window in which
% FILENAME does not exist is therefore one rename long (microseconds) instead of
% the several minutes a save() of a 200 MB model takes.
%
% CRASH MATRIX
% ------------
%   during the solve         final_md.mat untouched          re-run the leg
%   during save(tmp)         good final_md.mat + orphan tmp  G3 refuses the next
%                                                            run until the
%                                                            operator looks. THE
%                                                            SEED IS INTACT.
%   between steps 4 and 5    .leg<N>.mat present, final_md   the ONLY genuinely
%                            ABSENT                          bad window; recover
%                                                            with
%                                                            mv <stem>.leg<N>.mat
%                                                               <FILENAME>
%   after step 5             new file live, one generation   same recovery
%
% `relax.py --check <RUN>` detects and prints the recovery command for the two
% recoverable states.
%
% A GENERATION IS NOT A BACKUP. It survives one bad leg, not two, and it lives on
% the same Lustre filesystem in the same directory as the thing it protects. For
% an irreplaceable seed -- anything backing a shipped figure -- the right answer
% is a copy under a name no tool writes to. relax.py refuses in-place
% continuation on any run listed in Docs/protected_steady_states.txt unless
% --force is given, for exactly that reason.

    if nargin < 4 || isempty(keep_generations)
        keep_generations = 1;
    end
    if ~exist(run_dir, 'dir')
        error('ISSMMiscUtils.save_model_atomic: no such directory: %s', run_dir);
    end

    target = fullfile(run_dir, filename);
    [~, stem, ~] = fileparts(filename);

    % A job id makes the tmp name attributable: G3 refuses to start when one is
    % present, and the operator needs to know WHICH job left it behind.
    tag = getenv('SLURM_JOB_ID');
    if isempty(tag)
        tag = sprintf('pid%d', feature('getpid'));
    end
    % See the header: the trailing .mat is required for load() to work.
    tmp = fullfile(run_dir, sprintf('%s.tmp-%s.mat', stem, tag));

    if exist(tmp, 'file')
        error(['ISSMMiscUtils.save_model_atomic: %s already exists. Another ' ...
               'process with the same job id is mid-write, or a previous one ' ...
               'died there. Refusing to race it.'], tmp);
    end

    % ------------------------------------------------------------- step 2 ----
    % -v7.3 when the model is large: MAT-v5 caps a single variable at 2 GB and
    % save() fails at the END of a multi-minute write, i.e. in the worst place.
    % Sized from the in-memory model rather than guessed.
    w = whos('md');
    if w.bytes > 1.5e9
        fprintf(['[save_model_atomic] model is %.2f GB in memory -- writing ' ...
                 '-v7.3 (MAT-v5 caps one variable at 2 GB)\n'], w.bytes/1e9);
        save(tmp, 'md', '-v7.3');
    else
        save(tmp, 'md');
    end

    d = dir(tmp);
    if isempty(d) || d.bytes == 0
        error('ISSMMiscUtils.save_model_atomic: save() produced no bytes at %s', tmp);
    end

    % ------------------------------------------------------- step 3 (G14) ----
    % Reload and check the two things that would make the file useless. This is
    % the guard that catches a truncated or corrupt write BEFORE it replaces the
    % known-good file -- without it, "atomic rename" just means atomically
    % installing garbage.
    verified = false;
    try
        chk = load(tmp, 'md');
        rl  = chk.md;
        if rl.mesh.numberofvertices ~= md.mesh.numberofvertices
            error(['reloaded numberofvertices %d ~= %d in memory'], ...
                  rl.mesh.numberofvertices, md.mesh.numberofvertices);
        end
        if ~isfield(rl.results, 'TransientSolution') || isempty(rl.results.TransientSolution)
            error('reloaded model has no TransientSolution');
        end
        n_mem  = numel(md.results.TransientSolution);
        n_disk = numel(rl.results.TransientSolution);
        if n_mem ~= n_disk
            error('reloaded %d snapshots, %d in memory', n_disk, n_mem);
        end
        t_disk = double(rl.results.TransientSolution(end).time);
        t_mem  = double(md.results.TransientSolution(end).time);
        if abs(t_disk - t_mem) > 1e-6 * max(1, abs(t_mem))
            error('reloaded (end).time = %.6f, %.6f in memory', t_disk, t_mem);
        end
        if isempty(rl.results.TransientSolution(end).Thickness)
            error(['reloaded seed record has an EMPTY Thickness -- the seed ' ...
                   'must be a COMPLETE record (thin_transientsol never blanks ' ...
                   'it; something else did)']);
        end
        clear chk rl
        verified = true;
    catch ME
        % Leave the tmp on disk. It is evidence, and the good file is untouched.
        error(['ISSMMiscUtils.save_model_atomic: G14 VERIFY FAILED on %s (%s). ' ...
               'The existing %s has NOT been touched and is still good. The ' ...
               'bad tmp is left in place deliberately, as evidence; delete it ' ...
               'by hand once you have looked at it.'], tmp, ME.message, target);
    end

    % ------------------------------------------------------------- step 4 ----
    generation_path = '';
    if keep_generations > 0 && exist(target, 'file')
        next = next_generation_index(run_dir, stem);
        generation_path = fullfile(run_dir, sprintf('%s.leg%d.mat', stem, next));
        [ok, msg] = movefile(target, generation_path, 'f');
        if ~ok
            error(['ISSMMiscUtils.save_model_atomic: could not rotate %s -> %s ' ...
                   '(%s). Aborting BEFORE the swap; %s is still the good file ' ...
                   'and %s holds the new one.'], target, generation_path, msg, ...
                  target, tmp);
        end
    end

    % ------------------------------------------------------------- step 5 ----
    [ok, msg] = movefile(tmp, target, 'f');
    if ~ok
        % This is the bad window. Say exactly how to get out of it.
        error(['ISSMMiscUtils.save_model_atomic: FAILED to rename %s -> %s ' ...
               '(%s). YOU ARE IN THE BAD WINDOW: %s does not exist. Recover ' ...
               'with:\n    mv %s %s'], tmp, target, msg, target, ...
              generation_path, target);
    end

    % ------------------------------------------------------------- step 6 ----
    pruned = prune_generations(run_dir, stem, keep_generations);

    d = dir(target);
    info = struct('path', target, 'tmp_path', tmp, ...
                  'generation_path', generation_path, ...
                  'bytes', d.bytes, 'verified', verified, 'pruned', {pruned});
    fprintf('[save_model_atomic] wrote %s (%.1f MB, G14 verified)\n', ...
            target, d.bytes/1e6);
    if ~isempty(generation_path)
        fprintf('[save_model_atomic] previous generation kept at %s\n', ...
                generation_path);
    end
end


function next = next_generation_index(run_dir, stem)
% Generation indices are monotone and never reused, so a log line naming
% final_md.leg3.mat refers to one file for all time.
    listing = dir(fullfile(run_dir, sprintf('%s.leg*.mat', stem)));
    next = 1;
    for k = 1:numel(listing)
        tok = regexp(listing(k).name, '\.leg(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            next = max(next, str2double(tok{1}) + 1);
        end
    end
end


function pruned = prune_generations(run_dir, stem, keep_generations)
    pruned = {};
    listing = dir(fullfile(run_dir, sprintf('%s.leg*.mat', stem)));
    if isempty(listing); return; end
    idx = nan(1, numel(listing));
    for k = 1:numel(listing)
        tok = regexp(listing(k).name, '\.leg(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok); idx(k) = str2double(tok{1}); end
    end
    ok = ~isnan(idx);
    listing = listing(ok); idx = idx(ok);
    [~, order] = sort(idx);           % oldest generation first
    listing = listing(order);
    n_drop = numel(listing) - keep_generations;
    for k = 1:n_drop
        p = fullfile(run_dir, listing(k).name);
        delete(p);
        pruned{end+1} = p; %#ok<AGROW>
        fprintf('[save_model_atomic] pruned old generation %s\n', p);
    end
end
