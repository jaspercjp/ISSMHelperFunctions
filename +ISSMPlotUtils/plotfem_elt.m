function h = plotfem_elt(md, v, varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'View3D', false, @(x) islogical(x));
    addParameter(p, 'AMR', nan, @(x) isnan(x) || (isnumeric(x) && x>0));
    parse(p, varargin{:});
    
    VIEW_3D  = p.Results.View3D;
    AMR      = p.Results.AMR;
    Z_OFFSET = 500; 
    
    % 1. Extract Geometry
    if isa(md.mesh, 'mesh3dprisms')
        x=md.mesh.x2d; y=md.mesh.y2d; elements=md.mesh.elements2d;
        % Note: v should already be sized to md.mesh.numberofelements2d
    elseif AMR 
        x = md.results.TransientSolution(AMR).MeshX; 
        y = md.results.TransientSolution(AMR).MeshY;
        elements = md.results.TransientSolution(AMR).MeshElements;
    else
        x=md.mesh.x; y=md.mesh.y; elements=md.mesh.elements;
    end
    x=x/1e3; y=y/1e3;

    % 2. Validation
    num_elements = size(elements, 1);
    if length(v) ~= num_elements
        error('Input vector v (%d) must match the number of elements (%d).', length(v), num_elements);
    end

    % 3. Plotting
    % When v is [num_elements x 1], trisurf uses it for CData per face
    h = trisurf(elements, x, y, Z_OFFSET*ones(size(x)), v, ...
                'FaceColor', 'flat', ...
                'EdgeColor', 'none');

    % 4. View and Aesthetics
    if VIEW_3D
        view(-217, 21);
        camlight;
        lighting flat; % Use flat lighting for element-based data
    else
        view(2);
    end

    shading flat; % Critical for per-element visualization
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    colorbar; % Usually helpful for element data
end