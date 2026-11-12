function [cd_tw, H_tw, Rxx_tw, Rxx_IT_tw] = extract_width_int_cd(md, x_coords, t)

    % Get the floating vertices 
    % float_vert = md.results.TransientSolution(1).MaskOceanLevelset<=0.;
    % float_min_x = min(md.mesh.x(float_vert));
    % float_max_x = max(md.mesh.x(float_vert));

    % For each of the vertices, at each time step, calculate the 
    %     width-integrated stress and the crevasse depth from it 
    % N_COORDS = 100;
    % coords = linspace(float_min_x, float_max_x, N_COORDS);
    
    rho_i=md.materials.rho_ice; 
    rho_w=md.materials.rho_water; 
    g=md.constants.g;
    sxx_t = md.results.TransientSolution(t).DeviatoricStressxx; 
    syy_t = md.results.TransientSolution(t).DeviatoricStressyy; 
    H_t = md.results.TransientSolution(t).Thickness;
    Rxx_t = 2*sxx_t + syy_t; 
    % Rxx_IT_t = rho_i*g*(1-rho_i/rho_w)*H_t / 2.;   

    Rxx_tw = zeros(length(x_coords), 1);
    H_tw = zeros(length(x_coords), 1);
    Rxx_IT_tw = zeros(length(x_coords), 1);
    for i=1:length(x_coords)
        x = x_coords(i); 
        width_idx = abs(md.mesh.x - x)<=900; 
        y_transect = md.mesh.y(width_idx);
        [y_transect, sort_idx] = sort(y_transect);
        Rxx_transect=Rxx_t(width_idx); Rxx_transect=Rxx_transect(sort_idx);
        H_transect=H_t(width_idx); H_transect=H_transect(sort_idx);
        Rxx_tw(i) = trapz(y_transect, Rxx_transect) / (max(y_transect) - min(y_transect));
        H_tw(i) = trapz(y_transect, H_transect) / (max(y_transect) - min(y_transect));
        Rxx_IT_tw(i) = rho_i*g*(1-rho_i/rho_w)*H_tw(i)/2;
    end

    cd_tw = 1 - sqrt(max(1-Rxx_tw./Rxx_IT_tw, 0));
end