function animate_cd_ils(md, DATA_DIR, varargin)
    if nargin==3
        T_SKIP=varargin{1};
    else
        T_SKIP=5;
    end

    FRAME_DIR=DATA_DIR+"cd_ils_frames/";
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir %s", FRAME_DIR));
    end
    MIN_THICKNESS = -100;

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 21);
    cm_purp_green=cm_interp([0.7 1 0.35; 1 1 1; 1 0.2 0.6], 21);
    fig=figure; set(fig, 'Position', [100, 100, 700, 400]); % [left, bottom, width, height]
    tiledlayout(2,1,'TileSpacing','compact','Padding','loose');

    num_t = length(md.results.TransientSolution);
    for t=1:T_SKIP:num_t
        sol = md.results.TransientSolution(t);

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
        
        ils = sol.MaskIceLevelset; 
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
        % cd(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        plotfem(x/1000,y/1000,cd);
        clim([0 1.2]); colorbar; colormap(ax_tmp,jet(12));
        xticks([]); set(gca,'Color',[1 1 1]);
        % yticks([]); 
        title("Crevasse Depth"); ylabel('$y$ [km]');

        ax_tmp = nexttile(2);
        plotfem(x/1000,y/1000,ils);
        clim([-1 1]); colorbar; colormap(ax_tmp,jet(10));
        % xticks([]);yticks([]); 
        set(gca,'Color',[1 1 1]);
        title("Ice Levelset");
        xlabel('$x$ [km]'); ylabel('$y$ [km]');

        t_pert = sol.time - md.results.TransientSolution(1).time;
        
        sgt=sgtitle(sprintf("\\bf Time=%.1f yr", t_pert), ...
            'interpreter','latex', 'FontSize',20, 'Color', 'red');

        for i=[1 2]
            nexttile(i);
            xlim([500 700]);
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
                            FRAME_DIR, DATA_DIR + "cd_ils_animation.mp4", FRAME_RATE);
    system(command);
end