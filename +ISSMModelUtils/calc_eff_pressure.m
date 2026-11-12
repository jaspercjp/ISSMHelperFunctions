function [N, p_ice] = calc_eff_pressure(md)
%% computes the effective pressure assuming perfect ocean connection
    b = md.geometry.base;
    h = md.geometry.thickness;
    rho_ice = md.materials.rho_ice;
    rho_water = md.materials.rho_water;
    g = md.constants.g;
    p_ice = rho_ice * g * h;
    p_water = -rho_water * g * b;

    N = max(0, p_ice - p_water);
end
