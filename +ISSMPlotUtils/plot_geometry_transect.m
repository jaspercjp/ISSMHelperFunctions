function [slice_x, bed_slice, base_slice, surface_slice] = plot_geometry_transect(md, idx, y_trnsct)
%PLOT_GEOMETRY_TRANSECT Plot and return the bed / base / surface profiles along a
%   y = y_trnsct transect of `md` at TransientSolution index `idx`.
    assert(isfield(md.results, 'TransientSolution'), ...
           'The provided model must have a TransientSolution in its results field.');

    sol = md.results.TransientSolution(idx);
    base = sol.Base;
    surface = sol.Surface;

    bed = md.geometry.bed;

    [bed_slice, slice_x]     = ISSMModelUtils.extract_var_ytrnsct(md, bed, y_trnsct);
    [base_slice, slice_x]    = ISSMModelUtils.extract_var_ytrnsct(md, base, y_trnsct);
    [surface_slice, slice_x] = ISSMModelUtils.extract_var_ytrnsct(md, surface, y_trnsct);

    hold on;
    plot(slice_x, bed_slice, 'k-');
    plot(slice_x, surface_slice, 'b-');
    plot(slice_x, base_slice, 'b-');
end
