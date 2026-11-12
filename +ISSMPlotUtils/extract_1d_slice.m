function [slice_coord, slice_vals] = extract_1d_slice(x,y,vals,slice_loc,axis,varargin)
    if nargin==6
        threshold=varargin{1};
    else
        threshold=500;
    end

    % 1:x, 2:y
    if axis==1
        slice_nodes = abs(x - slice_loc) < threshold;
        slice_coord = y(slice_nodes);
    elseif axis==2
        slice_nodes = abs(y - slice_loc) < threshold;
        slice_coord = x(slice_nodes);
    else
        error('[extract_1d_slice]: Unsupported axis number')
    end
    slice_vals=vals(slice_nodes);
    [slice_coord, sort_idx]=sort(slice_coord);
    slice_vals=slice_vals(sort_idx);

end