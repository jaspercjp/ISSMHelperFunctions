function animate_H_stressratio(md, DATA_DIR ,varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'NumYrs', 200, @(x) isnumeric(x) && x>0);
    addParameter(p, 'Skip', 5, @(x) isnumeric(x) && x>0);
    addParameter(p, 'FrameRate', 20, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    NUM_YRS=p.Results.NumYrs
    SKIP=p.Results.Skip;
    FRAME_RATE=p.Results.FrameRate;

    FRAME_DIR=[DATA_DIR '/H_stressratio_frames/'];
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir %s", FRAME_DIR));
    end
    MIN_THICKNESS = 2;

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 51);
    cm_purp_green=cm_interp([0.7 1 0.35; 1 1 1; 1 0.2 0.6], 51);
    fig=figure; set(fig, 'Position', [100, 100, 900, 450]); % [left, bottom, width, height]
    tiledlayout(2,1,'TileSpacing','compact','Padding','loose');

    rho_i=md.materials.rho_ice; rho_w=md.materials.rho_water; 
    g=md.constants.g;

    num_t = length(md.results.TransientSolution);
    for t=1:SKIP:num_t
        sol = md.results.TransientSolution(t);
        t_pert = sol.time - md.results.TransientSolution(1).time;
        if (t_pert>NUM_YRS)
            break
        end
        
        vx1=sol.Vx; vy1=sol.Vy; H=sol.Thickness;
        sxx=sol.DeviatoricStressxx;
        syy=sol.DeviatoricStressyy;
        Rxx=2*sxx+syy;
        RxxIT=rho_i*g*(1-rho_i/rho_w)*H/2;
        vx1(H<MIN_THICKNESS)=nan;
        vy1(H<MIN_THICKNESS)=nan;
        sxx(H<MIN_THICKNESS)=nan;
        syy(H<MIN_THICKNESS)=nan;
        H(H<MIN_THICKNESS)=nan;

        if isempty(vx1) 
            continue; 
        end

        % H(sol.MaskIceLevelset>0)=nan;
        if isfield(sol, 'MeshX')
            AMR_t = t;
        else
            AMR_t = nan;
        end
        
        if isfield(sol, 'MaskIceLevelset')
            H(sol.MaskIceLevelset > 0)   = nan;
            Rxx(sol.MaskIceLevelset > 0) = nan;
        end
        ax_tmp=nexttile(1); ISSMPlotUtils.plotfem(md, H, 'AMR', AMR_t);
        clim([0 800]); colorbar; colormap(ax_tmp,parula);
        xticks([]);yticks([]); set(gca,'Color',[1 1 1]);
        title("$H$ [m]");

        ax_tmp=nexttile(2); ISSMPlotUtils.plotfem(md, Rxx./RxxIT, 'AMR', AMR_t);
        clim([0 1.2]); colorbar; colormap(ax_tmp,jet);
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
        elseif isa(md.basalforcings, 'basalforcingsismip6')
            floating_idx = find(md.results.TransientSolution(t).MaskOceanLevelset <= 0);
            melt_rate = mean(md.results.TransientSolution(t).BasalforcingsFloatingiceMeltingRate(floating_idx));
        elseif isa(md.basalforcings, 'basalforcingsismip6explicit') || isa(md.basalforcings, 'basalforcingsideal')
            floating_idx = find(md.results.TransientSolution(t).MaskOceanLevelset <= 0);
            if isfield(sol, 'BasalforcingsFloatingiceMeltingRate')
                melt_rate = median(md.results.TransientSolution(t).BasalforcingsFloatingiceMeltingRate(floating_idx));
            else
                melt_rate = 0;
            end
        else    
            melt_rate = md.basalforcings.floatingice_melting_rate(1);
        end

        sgt=sgtitle(sprintf("\\bf Time=%.1f yr, Melt Forcing=%.2f m/yr", t_pert, melt_rate), ...
            'interpreter','latex', 'FontSize',20, 'Color', 'red');

        for i=1:2
            nexttile(i); xlim([0, 320]);
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
                            FRAME_DIR, DATA_DIR + "H_stressratio_animation.mp4", FRAME_RATE);
    system(command);
end