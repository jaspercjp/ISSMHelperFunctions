function plot_transect_cd(md, x_loc, varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'Window', 1000, @(x) isnumeric(x) && x>0);
    addParameter(p, 'Skip', 5, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    WINDOW=p.Results.Window;
    SKIP=p.Results.Skip;

    rho_ice = md.materials.rho_ice; 
    rho_water = md.materials.rho_water;
    g=md.constants.g;
    delta=1-rho_ice/rho_water; 
    if isprop(md.basalforcings, 'floatingice_melting_rate')
        melt_rate = md.basalforcings.floatingice_melting_rate(1);
    elseif contains(class(md.basalforcings), 'linearbasalforcings')
        melt_rate = md.basalforcings.upperwater_melting_rate;
        deep_melt_rate = md.basalforcings.deepwater_melting_rate;
    elseif contains(class(md.basalforcings), 'basalforcingsismip6')
        melt_rate = nan;
    end
    acc_rate = md.smb.mass_balance(1);

    xq=x_loc; % vertical transect
    trnsct_idx = find(abs(md.mesh.x - xq)<=WINDOW);
    trnsct_x = md.mesh.x(trnsct_idx);
    trnsct_y = md.mesh.y(trnsct_idx);
    [trnsct_y,sort_idx] = sort(trnsct_y);
    trnsct_idx=trnsct_idx(sort_idx);
    syy_trnsct = {};
    sxx_trnsct = {};
    H_trnsct = {};
    time = [];
    syy = [];
    sxx = [];
    cnt = 0;
    MIN_THICKNESS=10;
    for t=1:SKIP:length(md.results.TransientSolution)
        cnt = cnt+1;
        disp(t)
        sol=md.results.TransientSolution(t); 
        time(cnt)=sol.time;
        H_tmp = sol.Thickness; H_tmp(sol.MaskIceLevelset>0 | sol.MaskOceanLevelset>0 | H_tmp<MIN_THICKNESS)=nan;
        sxx_tmp=sol.DeviatoricStressxx; sxx_tmp(sol.MaskIceLevelset>0 | sol.MaskOceanLevelset>0 | H_tmp<MIN_THICKNESS)=nan;
        syy_tmp=sol.DeviatoricStressyy; syy_tmp(sol.MaskIceLevelset>0 | sol.MaskOceanLevelset>0 | H_tmp<MIN_THICKNESS)=nan;
        cd_tmp=sol.CrevasseDepth;
        syy_trnsct{cnt} = syy_tmp(trnsct_idx);
        sxx_trnsct{cnt} = sxx_tmp(trnsct_idx);
        H_trnsct{cnt} = H_tmp(trnsct_idx);
        cd_trnsct{cnt} = cd_tmp(trnsct_idx);

        % if any(isnan(sxx_trnsct{cnt}) | isnan(syy_trnsct{cnt}) | isnan(H_trnsct{cnt}))
        %     break;
        % end
    end
    time=time-time(1);
    Rxx_norm_trnsct={};
    for i=1:cnt
        Rxx_trnsct_tmp = 2*sxx_trnsct{i} + syy_trnsct{i};
        Rxx_IT_trnsct_tmp = rho_ice*g*(1-rho_ice/rho_water)*H_trnsct{i}/2;
        Rxx_norm_trnsct{i} = Rxx_trnsct_tmp ./ Rxx_IT_trnsct_tmp;
        if all(isnan(Rxx_norm_trnsct{i})) 
            cnt=i;
            break;
        end
    end

    cmap=jet(cnt);
    disp(cmap)
    figure;
    tiledlayout(1,2,'padding','compact','tilespacing','compact');
    for i=1:cnt
        nexttile(1); hold on; grid on;
        cd_trnsct_matlab=1-sqrt(max(1-Rxx_norm_trnsct{i}, 0));
        plot(trnsct_y/1000, max(cd_trnsct_matlab, 0), 'Color', cmap(i,:), 'LineWidth', 1.5);

        nexttile(2); hold on; grid on;
        plot(trnsct_y/1000, max(cd_trnsct{i}, 0), 'Color', cmap(i,:), 'LineWidth', 1.5);
    end

    nexttile(2); ylim([0 1]); title('Crevasse Depth (ISSM Internal)')
    colormap(cmap); cb=colorbar;
    clim([time(1) time(end)]);  % Assuming 'i' indexes time from 1 to cnt
    ylabel(cb, 'Time [yr]', 'interpreter','latex');  % or actual time if available
    sgtitle(sprintf("\\bf Melt Rate=%.2fm/yr, Transect through $x$=%dkm", melt_rate, floor(xq/1000)), 'color','red', 'interpreter','latex');
    xlabel("$y$ [km]");

    nexttile(1); ylim([0 1]); title('Crevasse Depth (MATLAB Calculation)')
    ylabel("$\tilde d$");
    xlabel("$y$ [km]");

    set(gcf,'Position', [0, 0, 800, 400]);

    %% Plot Thickness and Equivalent Ice Tongue Stress
    figure; 
    tiledlayout(1,2,'padding','compact','tilespacing','compact');
    for i=1:cnt
        nexttile(1); hold on; grid on;
        plot(trnsct_y/1000, H_trnsct{i}, 'Color', cmap(i,:));

        nexttile(2); hold on; grid on;
        plot(trnsct_y/1000, rho_ice*g*(1-rho_ice/rho_water)*H_trnsct{i}/2/1e3, 'Color', cmap(i,:));
    end 
    sgtitle(sprintf("\\bf Melt Rate=%.2fm/yr, Transect through $x$=%dkm", melt_rate, floor(xq/1000)), 'color','red', 'interpreter','latex');
    
    nexttile(1); ylim([0 400]); title('Thickness [m]')
    ylabel('Thickness [m]');
    xlabel("$y$ [km]");

    nexttile(2); ylim([0, 200]); 
    ylabel("Equivalent Ice Tongue Stress [kPa]");
    xlabel("$y$ [km]");
    colormap(cmap); cb=colorbar;
    clim([time(1) time(end)]);  % Assuming 'i' indexes time from 1 to cnt
    ylabel(cb, 'Time [yr]', 'interpreter','latex');  % or actual time if available

    set(gcf,'Position', [800, 0, 800, 400]);
end