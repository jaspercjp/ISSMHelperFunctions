function animate_velocities_secPstress(md, DATA_DIR ,varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'NumYrs', 200, @(x) isnumeric(x) && x>0);
    addParameter(p, 'Skip', 5, @(x) isnumeric(x) && x>0);
    addParameter(p, 'FrameRate', 20, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    NUM_YRS=p.Results.NumYrs
    SKIP=p.Results.Skip;
    FRAME_RATE=p.Results.FrameRate;

    FRAME_DIR=[DATA_DIR '/velocity_secPstress_frames/'];
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir %s", FRAME_DIR));
    end
    MIN_THICKNESS = 35;

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 21);
    fig=figure; set(fig, 'Position', [100, 100, 900, 800]); % [left, bottom, width, height]
    tiledlayout(4,1,'TileSpacing','compact','Padding','loose');

    
    num_t = length(md.results.TransientSolution);
    for t=1:SKIP:num_t
        sol = md.results.TransientSolution(t);
        t_pert = sol.time - md.results.TransientSolution(1).time;
        if (t_pert>NUM_YRS)
            break
        end

        vx1=sol.Vx; vy1=sol.Vy; H=sol.Thickness;
        sxx=sol.DeviatoricStressxx; exx=sol.StrainRatexx; 
        syy=sol.DeviatoricStressyy; eyy=sol.StrainRateyy;
        sxy=sol.DeviatoricStressxy; exy=sxy./sxx.*exx;
        e_diff = exx-eyy; 
        s2 = (sxx+syy)/2 - sqrt((sxx-syy).^2/4 + sxy.^2);
        e_diff(H<MIN_THICKNESS)=nan;
        vx1(H<MIN_THICKNESS)=nan;
        vy1(H<MIN_THICKNESS)=nan;
        s2(H<MIN_THICKNESS) = nan;

        % vx1(sol.MaskIceLevelset>0)=nan;
        if isfield(sol, 'MeshX')
            AMR_t = t;
            if isempty(sol.MeshX)
                continue;
            end
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

        ax_tmp=nexttile(3); ISSMPlotUtils.plotfem(md,s2/1e3, 'AMR', AMR_t, 'ZOffset', 0);
        clim([-25 25]); colorbar; colormap(ax_tmp,cm_pert);
        title(" Second P. Dev. Stress [kPa]");

        ax_tmp=nexttile(4); hold on;
        yline(0, '-', 'LineWidth', 2, 'Color', [1 1 0]);
        ISSMPlotUtils.plotfem(md, e_diff, 'AMR', AMR_t);
        clim([-1 1]*1e-12); 
        colorbar; colormap(ax_tmp, cm_pert);
        title('$\dot\epsilon_{xx} - \dot\epsilon_{yy}$ [yr$^{-1}$]');

        if isa(md.basalforcings, 'basalforcingsismip6explicit') || isa(md.basalforcings, 'basalforcingsideal')
            floating_idx = find(md.results.TransientSolution(t).MaskOceanLevelset <= 0);
            if isfield(sol, 'BasalforcingsFloatingiceMeltingRate')
                melt_rate = median(md.results.TransientSolution(t).BasalforcingsFloatingiceMeltingRate(floating_idx));
            else
                melt_rate = 0;
            end
        else    
            melt_rate = md.basalforcings.floatingice_melting_rate(1);
        end

        sgt=sgtitle(sprintf("\\bf Time=%.1f yr, Avg. Basal Melt Rate=%.1f", t_pert, melt_rate), ...
            'interpreter','latex', 'FontSize',20, 'Color', 'red');

        % for i=1:4
        %     nexttile(i); xlim([0 320]);
        %     % xline(640, 'k--', 'LineWidth', 3)
        % end
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
                            FRAME_DIR, DATA_DIR + "velocity_animation.mp4", FRAME_RATE);
    system(command);
end