function Rxx_star = hydfrac_Rxx_star()
    K_Ic     = 150e3; % [Pa m^(1/2)]
    rho_i    = 917;
    g        = 9.8;
    f        = 1.068;
    F        = 1.122; 
    alpha    = (pi/3) * sqrt(2/3) * F^(3/2)/f^(1/2);
    Rxx_star = (K_Ic * sqrt(rho_i*g) / alpha)^(2/3); %[Pa]
        % QUESTION: is Rxx_star to be compared with depth-integrated Rxx or
        % surface Rxx? 
end