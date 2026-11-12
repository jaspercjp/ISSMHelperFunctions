% This script calculates the dimensionless eta and beta parameters that 
% quantify the relative scales of the viscous term to the lateral shear
% term
function [eta, beta] = calculate_eta_beta(md, steady_sol)
    B = md.materials.rheology_B(1); 
    rho=md.materials.rho_ice; rho_w=md.materials.rho_water; 
    g=md.constants.g; delta=(1-rho/rho_w);
    n=md.materials.rheology_n(1);

    if any(strcmp(fieldnames(steady_sol), 'MeshX'))
        x = steady_sol.MeshX;
        y = steady_sol.MeshY;
    else
        x = md.mesh.x;
        y = md.mesh.y; 
    end
    
    xx = linspace(min(x), max(x), 2000);
    yy = linspace(min(y), max(y), 2000);
    [xg, yg] = meshgrid(xx, yy);
    
    % interp spatial fields 
    ols_sip = scatteredInterpolant(x, y, steady_sol.MaskOceanLevelset);
        ols_interp = ols_sip(xg, yg);
    u_sip = scatteredInterpolant(x, y, steady_sol.Vx);
        u_interp = u_sip(xg, yg);
    h_sip = scatteredInterpolant(x, y, steady_sol.Thickness);
        h_interp = h_sip(xg, yg);

    gl_idx = find(abs(ols_interp)<5); 
        gl_y = yg(gl_idx); 
        [gl_y, sort_idx] = sort(gl_y); gl_idx=gl_idx(sort_idx); 
    ug = u_interp(gl_idx); 
    hg = h_interp(gl_idx);
    
    ug_wavg = mean(ug); ug_wavg=ug_wavg/md.constants.yts;
    hg_wavg = mean(hg);
    
    iceshelf_x = xg(ols_interp<=0);
    min_attach_x = min(xg(h_interp<5), [], 'all');
    if isempty(min_attach_x)
        min_attach_x=max(iceshelf_x, [], 'all');
    end
    Ls=min_attach_x - mean(xg(gl_idx), 'all');
    Ws=max(y) - min(y);

    qg_norm = hg_wavg * ug_wavg^(1/n) / (rho*delta*g*hg_wavg^2); 
    
    eta = 4*B*qg_norm / Ls^(1/n);
    beta = 2*(n+1)^(1/n) * B*Ls*qg_norm / Ws^(1/n+1);
    mdot_unit=Ls / (hg_wavg * ug_wavg * md.constants.yts);
    mdot=md.smb.mass_balance(1) - md.basalforcings.floatingice_melting_rate(1);
    
    fprintf('%11s = %3d found\n', 'GZ nodes', length(gl_idx));
    fprintf('%11s = %3d km\n', 'Ls', round(Ls/1000));
    fprintf('%11s = %3d km\n', 'Ws', round(Ws/1000));
    fprintf('%11s = %11.6f yr/m\n', 'mdot unit', mdot_unit);
    fprintf('%11s = %11.6f \n', 'nondim melt', mdot*mdot_unit);
    fprintf('%11s = %11.6f m/yr\n', 'ug', ug_wavg*md.constants.yts);
    fprintf('%11s = %11.6f m\n', 'hg', hg_wavg);
    fprintf('%11s = %11.6f\n', 'eta', eta);
    fprintf('%11s = %11.6f\n', 'beta', beta);
    fprintf('%11s = %11.6f\n', 'beta/eta', beta/eta);
end