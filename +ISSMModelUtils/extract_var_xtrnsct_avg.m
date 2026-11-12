% Assume that the y-direction is the longitudinal (a.k.a. along-flow) direction. 
function z_trnsct_avg = extract_var_xtrnsct_avg(md,var_field,x_slice,varargin)
    p = inputParser;
    addParameter(p, 'SearchL', 200, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    SEARCH_L=p.Results.SearchL;

    [z_slice, slice_y] = ISSMModelUtils.extract_var_xtrnsct(md, var_field, x_slice, 'SearchL', SEARCH_L);
    z_trnsct_avg = trapz(slice_y, z_slice) / (max(slice_y) - min(slice_y)); 
end