% Assume that the y-direction is the longitudinal (a.k.a. along-flow) direction. 
function [z_slice, x_slice_y] = extract_var_xtrnsct(md,var_field,x_slice,varargin)
    p = inputParser;
    addParameter(p, 'SearchL', 200, @(x) isnumeric(x) && x>0);
    addParameter(p, 'x', [], @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'y', [], @(x) isnumeric(x) || isempty(x));
    parse(p, varargin{:});
    SEARCH_L=p.Results.SearchL;
    X=p.Results.x;
    Y=p.Results.y;

    if ~isempty(X) && ~isempty(Y)
        x = X;
        y = Y;
    else
        x = md.mesh.x;
        y = md.mesh.y; 
    end

    x_slice_idx = find(abs(x-x_slice)<SEARCH_L);
    x_slice_y = y(x_slice_idx);
    [x_slice_y, sort_idx] = sort(x_slice_y);
    x_slice_idx=x_slice_idx(sort_idx);

    var_tmp = var_field(x_slice_idx); 
    z_slice=var_tmp;
end