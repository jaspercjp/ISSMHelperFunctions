function animate_velocities(md, DATA_DIR ,varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'NumYrs', 200, @(x) isnumeric(x) && x>0);
    addParameter(p, 'Skip', 5, @(x) isnumeric(x) && x>0);
    addParameter(p, 'FrameRate', 20, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    NUM_YRS=p.Results.NumYrs
    SKIP=p.Results.Skip;
    FRAME_RATE=p.Results.FrameRate;

    FRAME_DIR=DATA_DIR+"velocity_frames/";
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir %s", FRAME_DIR));
    end
    MIN_THICKNESS = 35;

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 21);
    fig=figure; set(fig, 'Position', [100, 100, 900, 400]); % [left, bottom, width, height]
    tiledlayout(2,1,'TileSpacing','compact','Padding','loose');

    
    num_t = length(md.results.TransientSolution);
    for t=1:SKIP:num_t
        sol = md.results.TransientSolution(t);
        t_pert = sol.time - md.results.TransientSolution(1).time;
        if (t_pert>NUM_YRS)
            break
        end

        if isfield(sol, 'MeshX')
            x=sol.MeshX; y=sol.MeshY;
        else 
            x=md.mesh.x; y=md.mesh.y; 
        end

        vx1=sol.Vx; vy1=sol.Vy; H=sol.Thickness;
        vx1(H<MIN_THICKNESS)=nan;
        vy1(H<MIN_THICKNESS)=nan;

        % if isempty(vx1) 
        %     continue; 
        % end

        % vx1(sol.MaskIceLevelset>0)=nan;
        if isfield(sol, 'MeshX')
            AMR_t = t;
        else
            AMR_t = nan;
        end
        ax_tmp=nexttile(1); ISSMPlotUtils.plotfem(md,vx1, 'AMR', AMR_t);
        clim([0 850]); colorbar; colormap(ax_tmp,parula(11));
        title("Along-Flow Velocity $u$ [m/yr]");
        
        % vy1(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(2); ISSMPlotUtils.plotfem(md,vy1, 'AMR', AMR_t);
        clim([-25 25]); colorbar; colormap(ax_tmp,cm_pert);
        title("Across-Flow Velocity $v$ [m/yr]");

        % if exist("config", 'var') && isfield(config, 'MELT_RATE_TIMES')
        %     if t_pert>=config.MELT_RATE_TIMES(1) && t_pert<config.MELT_RATE_TIMES(2)
        %         melt_rate = config.MELT_RATE(1);
        %     elseif t_pert>=config.MELT_RATE_TIMES(2) && t_pert<config.MELT_RATE_TIMES(3)
        %         melt_rate = config.MELT_RATE(2);
        %     elseif t_pert>=config.MELT_RATE_TIMES(3) && t_pert<config.MELT_RATE_TIMES(4)
        %         melt_rate = config.MELT_RATE(3);
        %     end
        % else
        %     melt_rate = md.basalforcings.floatingice_melting_rate(1);
        % end

        sgt=sgtitle(sprintf("\\bf Time=%.1f yr", t_pert), ...
            'interpreter','latex', 'FontSize',20, 'Color', 'red');

        for i=1:2
            nexttile(i); xlim([0 320]);
            % xline(640, 'k--', 'LineWidth', 3)
        end
         % Write to video 
        set(gcf, 'Color', 'w');            % Set figure background to white
        frame = getframe(gcf); % Capture current figure as a frame
        
        drawnow; 
        frame = getframe(gcf); % Capture current figure as a frame
        if t<100
            if t<10 
                frame_num="00"+string(t);
            else
                frame_num = "0" + string(t);
            end
        else
            frame_num = string(t);
        end
        imwrite(frame.cdata, FRAME_DIR+sprintf("frame%s.png",frame_num));  % Save as PNG image
    end
    command = sprintf("python3 %s %s %s %d", ISSMPlotUtils.write_to_video_script(), ...
                            FRAME_DIR, DATA_DIR + "velocity_animation.mp4", FRAME_RATE);
    system(command);
end