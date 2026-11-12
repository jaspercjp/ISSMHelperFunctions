function save_thin_amr_md(save_path, md)
%SAVE_THIN_AMR_MD  Write a resumability-only amr_md.mat (SAVE_AMR_MD = 'thin').
%
%   ISSMMiscUtils.save_thin_amr_md(SAVE_PATH, MD)
%
% MD is taken BY VALUE. MATLAB's copy-on-write means the thinning below cannot
% touch the caller's model -- which is load-bearing: calc_steady_state.m needs
% the FULL TransientSolution array immediately afterwards, to hand to
% ISSMMeshUtils.copy_amr_mesh.
%
% Saved under the variable name `md` so that every existing reader's
% `load(path, 'md')` keeps working unchanged. The file shrinks from 400-800 MB to
% ~10-20 MB; what it can still do (restore the refined mesh and the end state, so
% a crashed fix-icefront relax need not re-run the AMR transient) is unchanged.
% See ISSMMiscUtils.thin_amr_snapshots for why two records suffice.

    md.results.TransientSolution = ...
        ISSMMiscUtils.thin_amr_snapshots(md.results.TransientSolution);
    target = fullfile(save_path, 'amr_md.mat');
    save(target, 'md');
    d = dir(target);
    fprintf('[save_thin_amr_md] wrote %s (%.1f MB, %d snapshot(s))\n', ...
            target, d.bytes/1e6, numel(md.results.TransientSolution));
end
