function [sxx, sxy, syy, H, x_slice_y] = extract_stresses_xtrnsct(md,x_slice,t)
    % extract constants from the model
    rho_ice=md.materials.rho_ice; 
    rho_water=md.materials.rho_water;
    g=md.constants.g;
    MIN_THICKNESS=-2;

    x_slice_idx = find(abs(md.mesh.x-x_slice)<900);
    x_slice_y = md.mesh.y(x_slice_idx);
    [x_slice_y, sort_idx] = sort(x_slice_y);
    x_slice_idx=x_slice_idx(sort_idx);

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
    
    sxx = sxx_tmp(x_slice_idx);
    sxy = sxy_tmp(x_slice_idx);
    syy = syy_tmp(x_slice_idx);
    H = H_tmp(x_slice_idx);
    t = sol.time;
end