function view_steady_state(md)
   steady_sol = md.results.TransientSolution(1);
   if any(strcmp(fieldnames(steady_sol), 'MeshX'))
        x = steady_sol.MeshX; 
        y=steady_sol.MeshY; 
        elts = steady_sol.MeshElements;
    else
        x=md.mesh.x; 
        y=md.mesh.y;
        elts = md.mesh.elements;
    end

    figure; tiledlayout(3,1,'tilespacing','compact','padding','compact');
    nexttile; hold on;
    p = patch('Faces', elts, ...
          'Vertices', [x(:)/1000, y(:)/1000], ...
          'FaceColor', 'white', ...   % Color of the element body
          'EdgeColor', 'black', ...  % Color of the mesh lines
          'LineWidth', 0.5);

    nexttile; 
    plotfem(x/1e3,y/1e3,steady_sol.Vx); colorbar;

    ax_tmp=nexttile; ols_cmap=flipud(slanCM(53, 500));
    plotfem(x/1e3, y/1e3, steady_sol.MaskOceanLevelset); colorbar; clim([-1 1]);
    colormap(ax_tmp, ols_cmap);

    set(gcf,'Position', [100, 100, 600, 400]);

end