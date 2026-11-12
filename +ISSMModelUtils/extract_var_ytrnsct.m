% Assume that the y-direction is the longitudinal (a.k.a. along-flow) direction. 
function [var_slice, y_slice_x] = extract_var_ytrnsct(md,var_field,y_slice,varargin)
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
    
    % extract fields from the model
    y_slice_idx = find(abs(y-y_slice)<SEARCH_L);
    y_slice_x = x(y_slice_idx);
    [y_slice_x, sort_idx] = sort(y_slice_x);
    y_slice_idx=y_slice_idx(sort_idx);

    var_tmp = var_field;
    
    var_slice=var_tmp(y_slice_idx);
end