function contiguous_nodes = find_contiguous_calving_nodes(md, mask)
    % FIND_CONTIGUOUS_CALVING_NODES Find all nodes in one contiguous area adjacent to calving front
    %
    % Inputs:
    %   md   - ISSM model structure
    %   mask - Binary mask (md.mesh.numberofvertices x 1) where 1 = nodes of interest, 0 = other nodes
    %
    % Output:
    %   contiguous_nodes - Binary mask (md.mesh.numberofvertices x 1) where 1 = nodes in largest 
    %                     contiguous region adjacent to calving front, 0 = other nodes
    %
    % The calving front is defined as the zero levelset of md.mask.ice_levelset
    
    % Input validation
    if length(mask) ~= md.mesh.numberofvertices
        error('Mask size must be (md.mesh.numberofvertices, 1)');
    end
    
    % Ensure mask is binary
    mask = (mask > 0);
    
    % Find calving front nodes (ice levelset ~ 0)
    ice_levelset = md.mask.ice_levelset;
    levelset_tolerance = 1e-3;  % Adjust based on your mesh resolution
    
    % Calving front nodes: ice inside (levelset < 0) with neighbors outside (levelset > 0)
    calving_front_nodes = false(md.mesh.numberofvertices, 1);
    
    % Build adjacency list for efficient neighbor lookup
    adjacency = build_node_adjacency(md);
    
    % Find nodes that are at the ice-ocean boundary
    for i = 1:md.mesh.numberofvertices
        if ice_levelset(i) < levelset_tolerance  % This node is ice
            % Check if any neighbor is ocean (levelset > 0)
            neighbors = adjacency{i};
            if any(ice_levelset(neighbors) > -levelset_tolerance)
                calving_front_nodes(i) = true;
            end
        end
    end
    
    fprintf('Found %d calving front nodes\n', sum(calving_front_nodes));
    
    % Find masked nodes that are adjacent to calving front
    seed_nodes = false(md.mesh.numberofvertices, 1);
    
    for i = 1:md.mesh.numberofvertices
        if mask(i) && any(calving_front_nodes(adjacency{i}))
            seed_nodes(i) = true;
        end
    end
    
    if sum(seed_nodes) == 0
        warning('No masked nodes found adjacent to calving front');
        contiguous_nodes = false(md.mesh.numberofvertices, 1);
        return;
    end
    
    fprintf('Found %d seed nodes adjacent to calving front\n', sum(seed_nodes));
    
    % Find all connected components in the masked region
    visited = false(md.mesh.numberofvertices, 1);
    components = {};
    component_sizes = [];
    
    for i = 1:md.mesh.numberofvertices
        if mask(i) && ~visited(i)
            % Start BFS/DFS from this node
            component = flood_fill(i, mask, adjacency, visited);
            components{end+1} = component;
            component_sizes(end+1) = length(component);
            
            % Mark all nodes in this component as visited
            visited(component) = true;
        end
    end
    
    fprintf('Found %d connected components in mask\n', length(components));
    
    % Find which components contain seed nodes (adjacent to calving front)
    calving_adjacent_components = [];
    calving_component_sizes = [];
    
    for i = 1:length(components)
        component = components{i};
        if any(seed_nodes(component))
            calving_adjacent_components(end+1) = i;
            calving_component_sizes(end+1) = component_sizes(i);
        end
    end
    
    if isempty(calving_adjacent_components)
        warning('No connected components found adjacent to calving front');
        contiguous_nodes = false(md.mesh.numberofvertices, 1);
        return;
    end
    
    fprintf('Found %d components adjacent to calving front\n', length(calving_adjacent_components));
    
    % Return the largest component adjacent to calving front
    [~, max_idx] = max(calving_component_sizes);
    largest_component_idx = calving_adjacent_components(max_idx);
    largest_component = components{largest_component_idx};
    
    fprintf('Largest component has %d nodes\n', length(largest_component));
    
    % Create output mask
    contiguous_nodes = false(md.mesh.numberofvertices, 1);
    contiguous_nodes(largest_component) = true;
end

function adjacency = build_node_adjacency(md)
    % Build adjacency list for mesh nodes
    adjacency = cell(md.mesh.numberofvertices, 1);
    for i = 1:md.mesh.numberofvertices
        adjacency{i} = [];
    end
    
    % Add connections from elements
    elements = md.mesh.elements;
    for i = 1:size(elements, 1)
        elem_nodes = elements(i, :);
        elem_nodes = elem_nodes(elem_nodes > 0);  % Remove any padding zeros
        
        % Connect each node to all others in the element
        for j = 1:length(elem_nodes)
            for k = 1:length(elem_nodes)
                if j ~= k
                    node1 = elem_nodes(j);
                    node2 = elem_nodes(k);
                    adjacency{node1} = [adjacency{node1}, node2];
                end
            end
        end
    end
    
    % Remove duplicates and sort
    for i = 1:md.mesh.numberofvertices
        adjacency{i} = unique(adjacency{i});
    end
end

function component = flood_fill(start_node, mask, adjacency, visited)
    % Breadth-first search to find connected component
    component = [];
    queue = start_node;
    visited(start_node) = true;
    
    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];  % Remove first element
        component(end+1) = current;
        
        % Add unvisited masked neighbors to queue
        neighbors = adjacency{current};
        for neighbor = neighbors
            if mask(neighbor) && ~visited(neighbor)
                visited(neighbor) = true;
                queue(end+1) = neighbor;
            end
        end
    end
end