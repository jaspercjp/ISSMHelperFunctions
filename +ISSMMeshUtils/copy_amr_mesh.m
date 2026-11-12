function md = copy_amr_mesh(md, end_sol)
    % END_SOL may be a single TransientSolution snapshot or the whole array.
    %
    % ISSM writes the MeshX/MeshY/MeshElements triple onto a snapshot ONLY when
    % that save point coincides with an actual AMR remesh event. Whenever
    % md.settings.output_frequency and md.transient.amr_frequency are not
    % aligned, the LAST snapshot therefore frequently carries no mesh at all,
    % and `end_sol.MeshX` fails with "Unrecognized field name" -- the crash seen
    % on job 36121382. A uniform AMR metric (hmin == hmax, as used by the 1 km
    % configs) makes this the common case rather than the exception, because
    % bamg may decline to remesh at all.
    %
    % So: accept the whole array and scan BACKWARDS for the most recent snapshot
    % that actually carries a mesh. That snapshot holds the mesh the solve ended
    % on, which is what every caller wants.
    sol = [];
    for k = numel(end_sol):-1:1
        if isfield(end_sol(k), 'MeshX') && ~isempty(end_sol(k).MeshX)
            sol = end_sol(k);
            break
        end
    end
    if isempty(sol)
        % Deliberately fatal rather than a silent fallback: by the time this is
        % called the caller has usually run ReinitializeAMRModel, which rebuilds
        % md.mesh from the domain outline. Quietly "keeping the current mesh"
        % would hand back the COARSE mesh under the name of a refined one, which
        % is precisely the kind of silent resolution downgrade the steady-state
        % campaign exists to catch.
        error(['ISSMMeshUtils.copy_amr_mesh: none of the %d supplied ' ...
               'TransientSolution snapshots carries a MeshX/MeshY/MeshElements ' ...
               'triple, i.e. AMR never recorded a remesh in this solve. Check ' ...
               'that md.transient.amr_frequency > 0 and that the AMR metric ' ...
               '(md.amr.hmin/hmax/*_resolution) actually differs from the ' ...
               'current mesh.'], numel(end_sol));
    end

    % Update stats
    md.mesh.numberofvertices = numel(sol.MeshX);
    md.mesh.elements = sol.MeshElements;
    md.mesh.numberofelements = size(sol.MeshElements, 1);
    md.mesh.x = sol.MeshX;
    md.mesh.y = sol.MeshY;
    
    % Boundary vertices
    tol = 1e-2;
    md.mesh.vertexonboundary = abs(md.mesh.x-min(md.mesh.x))<tol | abs(md.mesh.x-max(md.mesh.x))<tol | ...
                               abs(md.mesh.y-min(md.mesh.y))<tol | abs(md.mesh.y-max(md.mesh.y))<tol;
    
    % Edges construction
    all_edges = [md.mesh.elements(:,[1,2]); md.mesh.elements(:,[2,3]); md.mesh.elements(:,[3,1])];
    [~, ia, ic] = unique(sort(all_edges, 2), 'rows');
    el_ids = repmat((1:md.mesh.numberofelements)', 3, 1);
    
    e1 = el_ids(ia);
    e2 = accumarray(ic, el_ids) - e1;
    e2(e2==0) = -1; % Boundary marker
    
    md.mesh.edges = [all_edges(ia,:), e1, e2];
    md.mesh.numberofedges = size(md.mesh.edges, 1);
    
    % Segments (boundary edges)
    is_seg = e2 == -1;
    md.mesh.segments = md.mesh.edges(is_seg, 1:3);
    md.mesh.segmentmarkers = ones(sum(is_seg), 1);
    
    % Element connectivity
    tr = triangulation(md.mesh.elements, md.mesh.x, md.mesh.y);
    md.mesh.elementconnectivity = neighbors(tr);
    md.mesh.elementconnectivity(isnan(md.mesh.elementconnectivity)) = 0;
    
    % Vertex connectivity
    md.mesh.average_vertex_connectivity = 9; 
end


function [edges, segments, segmentmarkers, elementconnectivity, avg_vertex_conn] = get_mesh_properties_issm(elements, num_vertices)
    % elements: Nx3 array of vertex indices
    % num_vertices: Total number of vertices (optional, inferred if not provided)
    
    if nargin < 2
        num_vertices = max(elements(:));
    end
    num_elements = size(elements, 1);
    
    %% 1. Calculate Edges
    % Create list of all edges from all elements
    e1 = elements(:, [1 2]);
    e2 = elements(:, [2 3]);
    e3 = elements(:, [3 1]);
    
    all_edges = [e1; e2; e3];
    all_elems = [1:num_elements, 1:num_elements, 1:num_elements]';
    
    % Sort vertices to identify undirected edges
    all_edges_sorted = sort(all_edges, 2);
    
    % Find unique edges
    [~, ia, ic] = unique(all_edges_sorted, 'rows');
    
    % Determine element1 (first element encountered)
    element1 = all_elems(ia);
    
    % Determine element2 using sum strategy
    sum_elems = accumarray(ic, all_elems);
    element2 = sum_elems - element1;
    
    % Set boundary edges to -1
    element2(element2 == 0) = -1;
    
    % Retrieve ordered vertices (CCW based on element1)
    v1 = all_edges(ia, 1);
    v2 = all_edges(ia, 2);
    
    edges = [v1, v2, element1, element2];
    
    %% 2. Calculate Segments and Markers
    % Segments are edges where element2 is -1 (boundary)
    boundary_mask = (element2 == -1);
    
    % Segments: [v1, v2, element1]
    segments = edges(boundary_mask, 1:3);
    
    % SegmentMarkers: Default to 1
    segmentmarkers = ones(size(segments, 1), 1);
    
    %% 3. Calculate Element Connectivity
    % Identify neighbors from internal edges (where element2 != -1)
    internal_mask = ~boundary_mask;
    el_pairs = edges(internal_mask, 3:4);
    
    % Create adjacency list (el1 <-> el2)
    % Stack both directions: [el1, el2] and [el2, el1]
    adjacency = [el_pairs; el_pairs(:, 2), el_pairs(:, 1)];
    
    % Sort by source element to group neighbors
    [sorted_src, sort_idx] = sort(adjacency(:, 1));
    sorted_neighbors = adjacency(sort_idx, 2);
    
    % Calculate column index (1, 2, or 3) for each neighbor
    % Finds transitions in sorted_src
    change_points = [1; diff(sorted_src) ~= 0];
    % Generate group IDs (1 to num_groups)
    group_ids = cumsum(change_points);
    % Generate 1-based index within each group
    % We compute 'start index' of each group and subtract from current index
    group_starts = accumarray(group_ids, (1:length(sorted_src))', [], @min);
    col_idx = (1:length(sorted_src))' - group_starts(group_ids) + 1;
    
    % Fill connectivity matrix
    elementconnectivity = zeros(num_elements, 3);
    lin_idx = sub2ind([num_elements, 3], sorted_src, col_idx);
    elementconnectivity(lin_idx) = sorted_neighbors;
    
    %% 4. Calculate Average Vertex Connectivity
    % Degree of a vertex is the number of edges connected to it
    % Using unique edges from 'edges' array (account for v1 and v2)
    degrees = accumarray([edges(:, 1); edges(:, 2)], 1, [num_vertices, 1]);
    
    % Average connectivity (graph degree)
    avg_vertex_conn = mean(degrees);
end