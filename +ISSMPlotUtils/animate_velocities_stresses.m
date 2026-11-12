function animate_dynamics(md, DATA_DIR ,varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'NumYrs', 200, @(x) isnumeric(x) && x>0);
    addParameter(p, 'Skip', 5, @(x) isnumeric(x) && x>0);
    addParameter(p, 'FrameRate', 20, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    NUM_YRS=p.Results.NumYrs
    SKIP=p.Results.Skip;
    FRAME_RATE=p.Results.FrameRate;

    FRAME_DIR=DATA_DIR+"velocity_stress_frames/";
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir %s", FRAME_DIR));
    end
    MIN_THICKNESS = 2;

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 21);
    cm_purp_green=cm_interp([0.7 1 0.35; 1 1 1; 1 0.2 0.6], 21);
    fig=figure; set(fig, 'Position', [100, 100, 900, 450]); % [left, bottom, width, height]
    tiledlayout(3,2,'TileSpacing','compact','Padding','loose');

    rho_i=md.materials.rho_ice; rho_w=md.materials.rho_water; 
    g=md.constants.g;

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
        sxx=sol.DeviatoricStressxx;
        syy=sol.DeviatoricStressyy;
        vx1(H<MIN_THICKNESS)=nan;
        vy1(H<MIN_THICKNESS)=nan;
        sxx(H<MIN_THICKNESS)=nan;
        syy(H<MIN_THICKNESS)=nan;

        if isempty(vx1) 
            continue; 
        end

        vx1(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(1); plotfem(x,y,vx1);
        clim([0 300]); colorbar; colormap(ax_tmp,parula(11));
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$v_x$ [m/yr]");
        
        vy1(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(2); plotfem(x,y,vy1);
        clim([-25 25]); colorbar; colormap(ax_tmp,cm_pert);
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$v_y$ [m/yr]");

        sxx(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(3); plotfem(x,y,sxx/1000);
        clim([-10 100]); colorbar; colormap(ax_tmp,parula(20));
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$\tau_{xx}$ [kPa]");

        syy(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(4); plotfem(x,y,syy/1000);
        clim([-20 20]); colorbar; colormap(ax_tmp,cm_pert);
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$\tau_{yy}$ [kPa]");

        Rxx=2*sxx+syy;
        ax_tmp=nexttile(5); plotfem(x,y,Rxx/1000);
        clim([0 100]); colorbar; colormap(ax_tmp,parula(20));
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$R_{xx}$ [kPa]");

        RxxIT=rho_i*g*(1-rho_i/rho_w)*H/2;
        ax_tmp=nexttile(6); plotfem(x,y,Rxx./RxxIT);
        clim([0 1]); colorbar; colormap(ax_tmp,jet(20));
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$R_{xx}/R_{xx}^{IT}$ [unitless]");

        if exist("config", 'var') && isfield(config, 'MELT_RATE_TIMES')
            if t_pert>=config.MELT_RATE_TIMES(1) && t_pert<config.MELT_RATE_TIMES(2)
                melt_rate = config.MELT_RATE(1);
            elseif t_pert>=config.MELT_RATE_TIMES(2) && t_pert<config.MELT_RATE_TIMES(3)
                melt_rate = config.MELT_RATE(2);
            elseif t_pert>=config.MELT_RATE_TIMES(3) && t_pert<config.MELT_RATE_TIMES(4)
                melt_rate = config.MELT_RATE(3);
            end
        else
            melt_rate = md.basalforcings.floatingice_melting_rate(1);
        end

        sgt=sgtitle(sprintf("\\bf Time=%.1f yr, Melt Forcing=%.2f m/yr", t_pert, melt_rate), ...
            'interpreter','latex', 'FontSize',20, 'Color', 'red');

        for i=1:6
            nexttile(i); xlim([5, 6.4]*1e5);
            set(gca,'Color',[0.8 0.8 0.8]);
        end
         
        % Write to video 
        set(gcf, 'Color', 'w');            % Set figure background to white
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
                            FRAME_DIR, DATA_DIR + "velocity_stress_animation.mp4", FRAME_RATE);
    system(command);
end