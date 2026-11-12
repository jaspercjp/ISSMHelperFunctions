function animate_dynamics(md, DATA_DIR ,varargin)
    if nargin==3
        T_SKIP=varargin{1};
    else
        T_SKIP=5;
    end

    FRAME_DIR=DATA_DIR+"dynamic_frames/";
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir %s", FRAME_DIR));
    end
    MIN_THICKNESS = -100;

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 21);
    cm_purp_green=cm_interp([0.7 1 0.35; 1 1 1; 1 0.2 0.6], 21);
    fig=figure; set(fig, 'Position', [100, 100, 900, 400]); % [left, bottom, width, height]
    tiledlayout(2,3,'TileSpacing','compact','Padding','loose');

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

        vx1(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(1); plotfem(x,y,vx1);
        clim([200 600]); colorbar; colormap(ax_tmp,parula(11));
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$v_x$ [m/yr]");
        
        vy1(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(4); plotfem(x,y,vy1);
        clim([-25 25]); colorbar; colormap(ax_tmp,cm_pert);
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$v_y$ [m/yr]");

        % Buttressing/Stress Plot
        ax_tmp=nexttile(3); 
        % if isfield(sol, 'ButtressingK')
        %     stress = sol.ButtressingK;
        %     stress(sol.MaskIceLevelset>0 | sol.Thickness<MIN_THICKNESS)=nan;
        %     plotfem(x,y,stress);
        %     title("Buttressing $K=1-\frac{2R_{yy}}{\rho g' H}$");
        %     clim([-1 1]); colorbar; colormap(ax_tmp, jet(21));
        % else
        %     Ryy_norm = (2*sxx+syy) ./ (rho_ice*g*(1-rho_ice/rho_water)*H/2);
        %     if t==1
        %         Ryy_norm_init = Ryy_norm; 
        %     end
        %     plotfem(x,y,Ryy_norm);
        %     clim([0 1.0]); colorbar; colormap(ax_tmp,jet(20));
        %     title("Norm Res. Stress, $\tilde R_{yy}$");
        % end
        ols(sol.MaskIceLevelset>0)=nan;
        plotfem(x,y,ols); clim([-1 1]); 
        title("Ocean Levelset"); colorbar;
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);

        % Thickness plot
        H(sol.MaskIceLevelset>0)=nan;
        ax_tmp=nexttile(5); plotfem(x,y,H);
        clim([300 600]);
        colorbar; colormap(ax_tmp,parula(15));
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("Thickness [m]")
        
        % Crevasse Depth plot
        ax_tmp = nexttile(2);
        if isfield(sol,'CrevasseDepth')
            cd=sol.CrevasseDepth./sol.Thickness;
        else
            rho_ice=md.materials.rho_ice; rho_water=md.materials.rho_water; 
            g=md.constants.g; delta=1-rho_ice/rho_water;
            K=1-(2*sol.DeviatoricStressyy+sol.DeviatoricStressxx) ./ (rho_ice*g*delta*sol.Thickness/2);
            cd=1-sqrt(max(K,0));
        end
        cd(sol.MaskIceLevelset>0)=nan;
        plotfem(x,y,cd);
        clim([0 1]);
        colorbar; colormap(ax_tmp,jet(10));
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("Crevasse Depth");
        
        % Ice levelset plot
        ax_tmp=nexttile(6); plotfem(x,y,ils);
        clim([-1 1])
        colorbar; colormap(ax_tmp,jet);
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("Ice Levelset")
        
        for l=1:6
            nexttile(l); %ylim([450000 1100000]);
        end

        t_pert = sol.time - md.results.TransientSolution(1).time;
        
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
                            FRAME_DIR, DATA_DIR + "dynamic_animation.mp4", FRAME_RATE);
    system(command);
end