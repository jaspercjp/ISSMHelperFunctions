% Assume that the y-direction is the longitudinal (a.k.a. along-flow) direction. 
function [sxx, sxy, syy, H, y_slice_x] = extract_stresses_ytrnsct(md,y_slice,t)
    % extract constants from the model
    rho_ice=md.materials.rho_ice; 
    rho_water=md.materials.rho_water;
    g=md.constants.g;
    MIN_THICKNESS=-2;

    y_slice_idx = find(abs(md.mesh.y-y_slice)<900);
    y_slice_x = md.mesh.x(y_slice_idx);
    [y_slice_x, sort_idx] = sort(y_slice_x);
    y_slice_idx=y_slice_idx(sort_idx);

    sol = md.results.TransientSolution(t);
    H_tmp = sol.Thickness;
    syy_tmp = sol.DeviatoricStressyy;
    sxx_tmp = sol.DeviatoricStressxx;
    sxy_tmp = sol.DeviatoricStressxy;
    
    % Ryy_tmp = 2*syy_tmp + sxx_tmp;
    % Ryy_IT_tmp = rho_ice*g*(1-rho_ice/rho_water)*H_tmp/2;

    % mask out time series when vertex has no ice
    if isfield(sol, 'MaskIceLevelset')
        invalid_idx=sol.MaskIceLevelset>0;
        H_tmp(invalid_idx)=nan;
        syy_tmp(invalid_idx)=nan;
        sxx_tmp(invalid_idx)=nan;
        sxy_tmp(invalid_idx)=nan;
        Ryy_tmp(invalid_idx)=nan;
        Ryy_IT_tmp(invalid_idx)=nan;
    end
    
    sxx = sxx_tmp(y_slice_idx);
    sxy = sxy_tmp(y_slice_idx);
    syy = syy_tmp(y_slice_idx);
    H = H_tmp(y_slice_idx);
    t = sol.time;
end