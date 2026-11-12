function out = concat_transientsol(a, b)
%CONCAT_TRANSIENTSOL  Concatenate two TransientSolution struct arrays whose
%   field sets may differ, by taking the union and filling the gaps with [].
%
%   OUT = ISSMMiscUtils.concat_transientsol(A, B)   % A first (older), then B
%
% WHY THIS IS NOT JUST `[a, b]`
% -----------------------------
% MATLAB refuses `[a, b]` when the two struct arrays do not have the SAME field
% names. Two things routinely make that happen across a leg boundary:
%
%   1. Thinning. A previously-thinned history keeps every field NAME (that is
%      the whole point of ISSMMiscUtils.thin_transientsol's blank-don't-rmfield
%      choice), so this case is usually fine -- but only usually.
%   2. AMR bookkeeping. ISSM writes the MeshX/MeshY/MeshElements triple onto a
%      snapshot ONLY when that save point coincides with an actual remesh event.
%      So one array can carry those three names and the other not, purely as an
%      accident of where the remeshes fell. This is the same asymmetry that
%      crashed job 36121382 through copy_amr_mesh.
%   3. A leg that ran with a different md.transient.requested_outputs. That is
%      supposed to be impossible (relax_continue never touches requested_outputs
%      and records its sha256 per leg so a violation is detectable after the
%      fact), but "supposed to be impossible" is not a reason to crash 8 hours
%      into a solve at the concatenation step rather than at the guard.
%
% MATLAB also tolerates differing field ORDER on concatenation, but relying on
% that is needless; orderfields is applied so the result is deterministic.
%
% The filler is [] rather than NaN because [] is what thin_transientsol writes
% for a dropped field, so a "field absent" and a "field thinned away" look
% identical downstream -- one empty-check covers both.

    if isempty(a); out = b; return; end
    if isempty(b); out = a; return; end

    fa = fieldnames(a);
    fb = fieldnames(b);
    only_b = fb(~ismember(fb, fa));
    only_a = fa(~ismember(fa, fb));

    for k = 1:numel(only_b)
        [a.(only_b{k})] = deal([]);
    end
    for k = 1:numel(only_a)
        [b.(only_a{k})] = deal([]);
    end

    a = orderfields(a);
    b = orderfields(b);
    out = [a(:)', b(:)'];
end
