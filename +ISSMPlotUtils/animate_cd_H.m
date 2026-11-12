function animate_cd_H(md, DATA_DIR, varargin)
    if nargin==3
        T_SKIP=varargin{1};
    else
        T_SKIP=5;
    end

    FRAME_DIR=DATA_DIR+"cd_H_frames/";
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir -p %s", FRAME_DIR));
    end
    MIN_THICKNESS = -10;

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 21);
    cm_purp_green=cm_interp([0.7 1 0.35; 1 1 1; 1 0.2 0.6], 21);
    fig=figure; set(fig, 'Position', [100, 100, 700, 400]); % [left, bottom, width, height]
    tiledlayout(2,1,'TileSpacing','compact','Padding','loose');

    num_t = length(md.results.TransientSolution);
    for t=2:T_SKIP:num_t
        sol = md.results.TransientSolution(t);
        last_sol = md.results.TransientSolution(t-1);

        if isfield(sol, 'MeshX')
            x=sol.MeshX; y=sol.MeshY;
        else 
            x=md.mesh.x; y=md.mesh.y; 
        end
        
        vx1=sol.Vx; vy1=sol.Vy; 
        syy=sol.DeviatoricStressyy;
        sxx=sol.DeviatoricStressxx;
        % vx1(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        % vy1(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        % syy(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        % sxx(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        vel1=sqrt(vx1.^2 + vy1.^2); 
        H=sol.Thickness;
        % H(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        
        ils = last_sol.MaskIceLevelset; 
        H = sol.Thickness; H(ils>0)=nan;
        ols = sol.MaskOceanLevelset;

        if isempty(vx1) || isempty(vel1)
            continue; 
        end
        
        % Crevasse Depth plot
        ax_tmp = nexttile(1);
        if isfield(sol,'CrevasseDepth')
            cd=sol.CrevasseDepth;
            cd(ils>0)=nan;
        else
            K=1-(2*sol.DeviatoricStressyy+sol.DeviatoricStressxx) ./ (rho_ice*g*delta*sol.Thickness/2);
            cd=1-sqrt(max(K,0));
        end
        cd(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        plotfem(x/1000,y/1000,cd);
        clim([0.0 1.1]); colorbar; colormap(ax_tmp,jet(11));
        ax_tmp.XTickLabel={}; 
        set(gca,'Color',[1 1 1]);
        % yticks([]); 
        title("Crevasse Depth"); ylabel('$y$ [km]'); grid on;

        ax_tmp = nexttile(2);
        H(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        plotfem(x/1000,y/1000,H);
        clim([100 800]); colorbar; colormap(ax_tmp,jet(10));
        % xticks([]);yticks([]); 
        set(gca,'Color',[1 1 1]);
        title("Thickness [m]");
        xlabel('$x$ [km]'); ylabel('$y$ [km]');  grid on;

        t_pert = sol.time - md.results.TransientSolution(1).time + md.timestepping.time_step;
        
        sgt=sgtitle(sprintf("\\bf Time=%.1f yr", t_pert), ...
            'interpreter','latex', 'FontSize',20, 'Color', 'red');

        for i=[1 2]
            nexttile(i);
            xlim([480 650]);
            xline(640, 'b--', 'LineWidth', 2);
            % xline(500, 'r-', 'LineWidth', 4);
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
    FRAME_RATE = 40;
    command = sprintf("python3 %s %s %s %d", ISSMPlotUtils.write_to_video_script(), ...
                            FRAME_DIR, DATA_DIR + "cd_H_animation.mp4", FRAME_RATE);
    system(command);
end