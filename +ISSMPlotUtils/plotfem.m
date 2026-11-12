function h = plotfem(md, v, varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'Color', 'red', @(x) ischar(x));
    addParameter(p, 'View3D', false, @(x) islogical(x));
    addParameter(p, 'AMR', nan, @(x) isnan(x) || (isnumeric(x) && x>0));
    addParameter(p, 'ZOffset', 0, @(x) isnumeric(x));
    parse(p, varargin{:});
    COLOR    = p.Results.Color
    VIEW_3D  = p.Results.View3D;
    AMR      = p.Results.AMR;
    Z_OFFSET = p.Results.ZOffset; % Use when overlaying images 
    
    if isa(md.mesh, 'mesh3dprisms')
        x=md.mesh.x2d; y=md.mesh.y2d; elements=md.mesh.elements2d;

        % handle the case where the input is a three-dimensional field
        if length(v) > md.mesh.numberofvertices2d
            v = v(1 : md.mesh.numberofvertices2d);
        end
    elseif ~isnan(AMR) 
        x = md.results.TransientSolution(AMR).MeshX; 
        y = md.results.TransientSolution(AMR).MeshY;
        elements = md.results.TransientSolution(AMR).MeshElements;
    else
        x=md.mesh.x; y=md.mesh.y; elements=md.mesh.elements;
    end
    x=x/1e3; y=y/1e3;

    % 1. Plotting using the explicit connectivity matrix
    % 'elements' should be an [N x 3] matrix of vertex indices
    if VIEW_3D
        h = trisurf(elements, x, y, v-Z_OFFSET, v, 'FaceColor', COLOR, 'EdgeColor', 'none');
        view(-217, 21);
        lighting gouraud;
        material shiny;
    else
        h = trisurf(elements, x, y, v-Z_OFFSET, v, 'FaceColor', 'interp', 'EdgeColor', 'none');
        view(2);
    end
    
    shading interp;
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
end