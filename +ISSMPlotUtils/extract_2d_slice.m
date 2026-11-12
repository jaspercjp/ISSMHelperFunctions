function [slice_coord_1, slice_coord_2, slice_vals] = extract_2d_slice(md,sol_name,var,slice_loc,axis,varargin)
    if nargin==6
        threshold=varargin{1};
    else
        threshold=500;
    end
    vals = md.results.(sol_name).(var);
    
    % 1:x, 2:y, 3:z
    if axis==1
        slice_nodes = abs(md.mesh.x - slice_loc) < threshold;
        slice_coord_1 = md.mesh.y(slice_nodes);
        slice_coord_2 = md.mesh.z(slice_nodes);
    elseif axis==2
        slice_nodes = abs(md.mesh.y - slice_loc) < threshold;
        slice_coord_1 = md.mesh.x(slice_nodes);
        slice_coord_2 = md.mesh.z(slice_nodes);
    elseif axis==3 
        slice_nodes = abs(md.mesh.z - slice_loc) < threshold;
        slice_coord_1 = md.mesh.x(slice_nodes);
        slice_coord_2 = md.mesh.y(slice_nodes);
    end
    slice_vals=vals(slice_nodes);

end