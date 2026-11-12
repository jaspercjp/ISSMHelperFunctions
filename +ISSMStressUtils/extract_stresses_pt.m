% Assume that the y-direction is the longitudinal (a.k.a. along-flow) direction. 
function [Ryy, Ryy_IT, sxx, syy, H, t] = extract_stresses_pt(md,x_pt,y_pt)
    % extract constants from the model
    rho_ice=md.materials.rho_ice; 
    rho_water=md.materials.rho_water;
    g=md.constants.g;
    MIN_THICKNESS=-2;

    % first find a number of x coordinates close to desired x point 
    dist_to_x_pt = abs(md.mesh.x - x_pt);
    x_isvalid = abs(dist_to_x_pt-min(dist_to_x_pt))<600;

    % from those points, find the closest to the desired y point 
    candidate_idx = find(x_isvalid); min_y_idx=inf; min_y_dist=inf;
    candidate_idx
    for i=1:length(candidate_idx)
        idx=candidate_idx(i);
        y_dist = abs(md.mesh.y(idx) - y_pt);
        if y_dist<min_y_dist
            min_y_idx=idx
            min_y_dist=y_dist;
        end
    end
    vert_idx=min_y_idx

    % use the physical coordinate values to determine the vertex number 
    x = md.mesh.x(vert_idx);
    y = md.mesh.y(vert_idx);
    [x,y]
    
    % extract the stresses at the point 
    num_steps = length(md.results.TransientSolution);
    t = nan*ones(1,num_steps);
    H = nan*ones(1,num_steps);
    sxx=nan*ones(1,num_steps);
    syy=nan*ones(1,num_steps);
    Ryy = nan*ones(1,num_steps);
    Ryy_IT = nan*ones(1,num_steps);
    Ryy_confined=nan*ones(1,num_steps);
    Ryy_IT_confined=nan*ones(1,num_steps);
    sxx = nan*ones(1,num_steps);
    syy = nan*ones(1,num_steps);
    for i=1:num_steps
        sol = md.results.TransientSolution(i);
        H_tmp = sol.Thickness;
        syy_tmp = sol.DeviatoricStressyy;
        sxx_tmp = sol.DeviatoricStressxx;
        
        Ryy_tmp = 2*syy_tmp + sxx_tmp;
        Ryy_IT_tmp = rho_ice*g*(1-rho_ice/rho_water)*H_tmp/2;

        % mask out time series when vertex has no ice
        if isfield(sol, 'MaskIceLevelset')
            invalid_idx=sol.MaskIceLevelset>0 ;
            H_tmp(invalid_idx)=nan;
            syy_tmp(invalid_idx)=nan;
            sxx_tmp(invalid_idx)=nan;
            Ryy_tmp(invalid_idx)=nan;
            Ryy_IT_tmp(invalid_idx)=nan;
        end

        Ryy(i) = Ryy_tmp(vert_idx);
        Ryy_IT(i) = Ryy_IT_tmp(vert_idx);
        sxx(i) = sxx_tmp(vert_idx);
        syy(i) = syy_tmp(vert_idx);
        H(i) = H_tmp(vert_idx);
        sxx(i) = sxx_tmp(vert_idx); 
        syy(i) = syy_tmp(vert_idx);
        t(i) = sol.time;
    end
    t=t-t(1);
end