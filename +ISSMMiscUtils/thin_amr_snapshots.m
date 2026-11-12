function sols = thin_amr_snapshots(sols)
%THIN_AMR_SNAPSHOTS  Reduce an AMR-phase TransientSolution to the only two
%   snapshots that amr_md.mat is ever actually used for.
%
%   SOLS = ISSMMiscUtils.thin_amr_snapshots(SOLS)
%
% ---------------------------------------------------------------------------
% WHY amr_md.mat CAN BE THINNED TO TWO RECORDS
% ---------------------------------------------------------------------------
% At 400-800 MB it is the single largest artefact calc_steady_state writes -- 5.9
% GB across the ~20 run folders on disk -- and NOTHING ON THE LIVE PATH READS IT
% FROM DISK. A repo-wide search over Scripts/, Notebooks/ and Docs/ finds it only
% in (a) the save() that writes it, (b) PHASE_FILES in steady_state_gate.py,
% which explicitly REFUSES to grade the AMR phase (mean(H) jumps at every remesh
% and check_steady_state's SkipAMR is inoperative on every config here), and
% (c) two throwaway scratch scripts. The AMR mesh reaches the next phase IN
% MEMORY via ISSMMeshUtils.copy_amr_mesh, never re-read from disk, and
% compare_steady_states._load_phases loads only `coarse` and `final`.
%
% Its one real value is RESUMABILITY: if the fix-icefront relax crashes, having
% amr_md.mat means not re-running the whole AMR transient. That is worth keeping,
% and it costs exactly two records:
%
%   1. the LAST snapshot -- the state the phase ended on;
%   2. the last snapshot CARRYING a MeshX/MeshY/MeshElements triple -- which is
%      what copy_amr_mesh scans backwards for.
%
% Those are frequently the same record, in which case one is kept.
%
% WHY (2) IS NOT JUST (1). ISSM writes the mesh triple onto a snapshot ONLY when
% that save point coincides with an actual AMR remesh event. Whenever
% output_frequency and amr_frequency are unaligned the last snapshot carries no
% mesh at all -- and a UNIFORM AMR metric (the 1 km campaign) makes that the
% common case rather than the exception, because bamg may decline to remesh once
% the mesh already matches the metric. That asymmetry is exactly what crashed job
% 36121382. Keeping both records means a resumed run can call copy_amr_mesh on
% the saved array and get the same answer it would have got in memory.
%
% This function is CONSERVATIVE BY CONSTRUCTION: if no snapshot carries a mesh it
% keeps the last one and says so, rather than erroring -- unlike copy_amr_mesh,
% which must be fatal there because it is about to hand back a mesh. Here the
% worst case is a slightly less resumable archive of a file nothing reads.

    n = numel(sols);
    if n <= 2
        return
    end

    keep = n;                      % (1) the last snapshot
    mesh_idx = [];
    for k = n:-1:1
        if isfield(sols(k), 'MeshX') && ~isempty(sols(k).MeshX)
            mesh_idx = k;
            break
        end
    end
    if isempty(mesh_idx)
        fprintf(2, ['[thin_amr_snapshots] WARNING: none of the %d AMR ' ...
                    'snapshots carries a MeshX triple, so the thinned ' ...
                    'amr_md.mat cannot be used to restore the refined mesh. ' ...
                    'Keeping the last snapshot only.\n'], n);
    elseif mesh_idx ~= n
        keep = [mesh_idx, n];
    end

    fprintf('[thin_amr_snapshots] %d snapshots -> %d (indices %s)\n', ...
            n, numel(keep), mat2str(keep));
    sols = sols(keep);
end
