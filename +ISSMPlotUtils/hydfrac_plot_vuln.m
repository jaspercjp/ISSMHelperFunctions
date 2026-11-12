function hydfrac_plot_vuln(md, vulnerability_mask, varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'DisplayLegend', false, @(x) islogical(x));
    parse(p, varargin{:});
    DISPLAY_LEGEND = p.Results.DisplayLegend;

    BG_COLOR = [0.9 0.9 0.9];
    SAFE_COLOR = [105, 172, 240] / 255;
    VULN_COLOR = [240, 105, 105] / 255;

    vulnerability_mask(md.mask.ice_levelset >= 0) = nan;
    if isfield(md.mesh, 'x2d')
        h=ISSMPlotUtils.plotfem(md, vulnerability_mask);
    else
        h=ISSMPlotUtils.plotfem(md, vulnerability_mask);
    end
    % alphaMask = double(~isnan(vulnerability_mask)); 
    % alphaMask = alphaMask(:); 
    % set(h, 'FaceColor', 'interp', ...      % or 'flat' depending on your needs
    %        'FaceAlpha', 'interp', ...      % Enables vertex-level transparency
    %        'FaceVertexAlphaData', alphaMask, ... 
    %        'AlphaDataMapping', 'none');    % Interprets 0 as transparent, 1 as opaque
    
    set(gca, 'Color', BG_COLOR);
    colormap([SAFE_COLOR; VULN_COLOR]);  clim([-1 1]);
    % colormap jet(20); clim([-1 1]);
    xlabel('$x$ [km]'); ylabel('$y$ [km]');
    grid off;

    if DISPLAY_LEGEND
        h1 = line(NaN, NaN, 'LineWidth', 2, 'Color', 'none', ...
            'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', VULN_COLOR, 'MarkerEdgeColor', VULN_COLOR); % Red
        h2 = line(NaN, NaN, 'LineWidth', 2, 'Color', 'none', ...
            'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', SAFE_COLOR, 'MarkerEdgeColor', SAFE_COLOR); % Green
        h3 = line(NaN, NaN, 'LineWidth', 2, 'Color', 'none', ...
            'Marker', 's', 'MarkerSize', 10, 'MarkerFaceColor', BG_COLOR, 'MarkerEdgeColor', BG_COLOR); % Gray
    
        lgd=legend([h1, h2, h3], {'Vulnerable ($R_{xx}>R_{xx}^*$)', 'Stable', 'Grounded Ice / No Ice'}, ...
            'Location', 'northoutside', 'Box', 'on');
        set(lgd, 'Color', 'w', 'TextColor', 'k');
    end
end