function [gl_x, gl_y, qx_g, qy_g] = calc_gl_flux(x, y, sol, varargin)
    p = inputParser;
    addParameter(p, 'SearchL', 200, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    SEARCH_L=p.Results.SearchL;

    u=sol.Vx; v=sol.Vy; H=sol.Thickness;
    gl_idx = find(abs(sol.MaskOceanLevelset)<SEARCH_L);
    gl_y = y(gl_idx); gl_x=x(gl_idx);
    [gl_y, sort_idx] = sort(gl_y);
    gl_idx=gl_idx(sort_idx);
    u_g=u(gl_idx); v_g=v(gl_idx); h_g=H(gl_idx); 
    qx_g=u_g.*h_g; qy_g=v_g.*h_g;

end