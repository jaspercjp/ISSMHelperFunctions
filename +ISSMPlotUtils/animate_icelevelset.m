function animate_icelevelset(md, DATA_DIR)
    videoFile = DATA_DIR+sprintf("ice_levelset_animation.mp4"); % Output file name
    v = VideoWriter(videoFile, 'MPEG-4'); % Create video writer object (MP4 format)
    v.FrameRate = 50;
    open(v); % Open video file for writing
    % MIN_THICKNESS = 5;
    cm_pert = cm_interp([0 0 1; 1 1 1; 1 0 0], 500);
    fig=figure; set(fig, 'Position', [100, 100, 400, 500]); % [left, bottom, width, height]
    tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    set(gcf,'Position',[100 100 800 200]);

    num_t = length(md.results.TransientSolution);
    for t=1:num_t
        sol = md.results.TransientSolution(t);

        if isfield(sol, 'MeshX')
            x=sol.MeshX; y=sol.MeshY;
        else 
            x=md.mesh.x; y=md.mesh.y; 
        end

        ils1 = sol.MaskIceLevelset;

        ax_tmp=nexttile(1); plotfem(x,y,ils1); shading flat; colorbar;
        clim([-1 1]); colormap(ax_tmp,jet); set(gca,'Color',[1 1 1]);
        xticks([]);yticks([]); 
        title(sprintf("Ice Shelf Geometry %.1f Years", sol.time-md.results.TransientSolution(1).time-md.timestepping.time_step));

        % Write to video 
        set(gcf, 'Color', 'w');      % Set figure background to white
        frame = getframe(gcf);       % Capture current figure as a frame
        writeVideo(v, frame);        % Write frame to video
    end
    close(v);
    disp('Ice levelset animation saved.');

end