
% compute_floating_area tokes in an ISSM model containing results from a transient simulation
% and computes the floating area with and without grounding line migration
function [a_float, a_float_init_adj, t] = calc_floating_area(md)
    num_steps = length(md.results.TransientSolution);

    init_ocean_mask = md.results.TransientSolution(1).MaskOceanLevelset;
    init_ice_mask=md.results.TransientSolution(1).MaskIceLevelset;
    init_float_elt=ISSMAreaUtils.extract_floating_elements(md, init_ocean_mask, init_ice_mask);

    % setup output variables 
    t = nan*ones(1,num_steps);
    a_float = nan*ones(1,num_steps);
    a_float_init_adj=nan*ones(1,num_steps);

    for i=1:num_steps
        sol = md.results.TransientSolution(i);
        % check which vertices are now floating 
        now_float_vert=find(sol.MaskOceanLevelset<0 & sol.MaskIceLevelset<0);
        now_float_elt={}; count=0;
        for j=1:length(init_float_elt)
            % check which initially floating vertices are still floating
            if any(init_float_elt{j}(1)==now_float_vert) && any(init_float_elt{j}(2)==now_float_vert) && any(init_float_elt{j}(3)==now_float_vert)
                count=count+1;
                now_float_elt{count}=init_float_elt{j};
            end
        end
        a_float_init_adj(i)=ISSMAreaUtils.calc_area_from_elt(md, now_float_elt);
        a_float(i) = sol.FloatingArea;
        t(i) = sol.time;
    end
    t=t-t(1);
end