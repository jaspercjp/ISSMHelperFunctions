%CHECK_STEADY_STATE Checks for steady state in ISSM model results over a specified period.
%
%   [time, MS_diff_H] = CHECK_STEADY_STATE(md)
%   [time, MS_diff_H] = CHECK_STEADY_STATE(md, 'NumYrs', NUM_YEARS)
%   [time, MS_diff_H, pass] = CHECK_STEADY_STATE(md, 'Plot', false, 'Threshold', 1e-3)
%
%   This function analyzes the last NumYrs of an ISSM model simulation to determine
%   whether the ice thickness has reached a steady state. It computes the mean
%   thickness change between consecutive output steps and normalizes it by the mean
%   surface mass balance, so MS_diff_H is dimensionless: 1.0 means the mean
%   thickness is changing by one year's accumulation per year.
%
%   Inputs:
%       md          - ISSM model structure containing simulation results and settings.
%   Name-value pairs:
%       'NumYrs'    - Window (yr) used for the plot's xlim and for the pass test
%                     (default 300).
%       'Plot'      - Draw the figure (default TRUE, so every existing interactive
%                     call is unchanged). Pass FALSE to use this in a batch job:
%                     figure needs a display (the log-scale call itself is
%                     version-safe: set(gca,'YScale','log'), not yscale()).
%       'Threshold' - If given, the third output `pass` is TRUE when the terminal
%                     MS_diff_H is at or below it AND the series is not increasing
%                     over the last NumYrs. Empty (default) leaves `pass` empty.
%       'SkipAMR'   - Skip output steps that coincide with an AMR remesh
%                     (default false). See the note below -- on most configs in
%                     this repo the test is degenerate and is auto-disabled.
%       'CutoffX'   - INERT. Parsed for backward compatibility only; its single
%                     consumer inside the loop has been commented out for as long
%                     as this file has existed.
%
%   Outputs:
%       time      - Array of time values for the checked period.
%       MS_diff_H - Array of normalized mean thickness differences.
%       pass      - Logical verdict against 'Threshold', or [] if none was given.
%
%   The function generates ONE plot: normalized mean thickness change over time.
%   (An earlier version also plotted grounding-line flux; that subplot is
%   commented out, so no grounding-line quantity is analysed here despite what
%   this docstring used to claim. Use ISSMGLUtils.calc_gl_flux_t for that.)
%
%   NOTE ON 'SkipAMR': the skip test is
%   mod(t*md.settings.output_frequency, md.transient.amr_frequency)==0, which is
%   only meaningful when the two frequencies are incommensurate. Every config in
%   HFBCalvingParametrization has output_frequency an exact multiple of
%   amr_frequency (5 | 25, 5 | 125, 5 | 250), so the test would skip 100% of
%   output steps and return an EMPTY series. That case is now detected, warned
%   about and disabled rather than silently returning no data.
%
%   The headless equivalent of this function, which reproduces MS_diff_H exactly
%   from a saved .mat without MATLAB, is
%   Scripts/Workflows/SteadyState/Sbatch/steady_state_gate.py -- it is what the
%   validation campaign gates on. See Docs/ADRs/0015 for the calibrated bars.
%
%   Example:
%       [time, MS_diff_H] = CHECK_STEADY_STATE(md, 'NumYrs', 500)
%
%   See also: ISSMGLUtils.extract_gl_flux, ISSMGLUtils.calc_gl_flux_t

function [time, MS_diff_H, pass] = check_steady_state(md, varargin)
    disp('[check_steady_state] :: Validating steady state ...')
    p = inputParser;

    % model setup
    addParameter(p, 'NumYrs', 300, @(x) isnumeric(x) && x>0);
    addParameter(p, 'CutoffX', 0, @(x) isnumeric(x) && x>0);
    addParameter(p, 'SkipAMR', false, @(x) islogical(x));
    addParameter(p, 'Plot', true, @(x) islogical(x));
    addParameter(p, 'Threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
    parse(p, varargin{:});
    NUM_YEARS=p.Results.NumYrs;
    CUTOFF_X=p.Results.CutoffX; %#ok<NASGU> inert; kept for call compatibility
    SKIP_AMR=p.Results.SkipAMR;
    DO_PLOT=p.Results.Plot;
    THRESHOLD=p.Results.Threshold;
    pass=[];

    % The SkipAMR test skips every output step when output_frequency is a
    % multiple of amr_frequency, which returns an empty series with no warning.
    % Detect that and fall back to keeping every step.
    if SKIP_AMR
        if md.transient.amr_frequency<=0
            warning('check_steady_state:SkipAMRUnused', ...
                ['SkipAMR requested but md.transient.amr_frequency = %g; ' ...
                 'no step can be skipped.'], md.transient.amr_frequency);
            SKIP_AMR=false;
        elseif mod(md.settings.output_frequency, md.transient.amr_frequency)==0
            warning('check_steady_state:SkipAMRDegenerate', ...
                ['SkipAMR disabled: output_frequency (%g) is a multiple of ' ...
                 'amr_frequency (%g), so the skip test would drop 100%% of ' ...
                 'output steps and return an empty series.'], ...
                md.settings.output_frequency, md.transient.amr_frequency);
            SKIP_AMR=false;
        end
    end

    num_steps=round(NUM_YEARS/(md.timestepping.time_step*md.settings.output_frequency));
    length_sim=length(md.results.TransientSolution);

    % check thickness 
    MS_diff_H=nan*ones(length_sim, 1);
    mean_H=nan*ones(length_sim,1);
    mean_Hg=nan*ones(length_sim,1);
    mean_ug=nan*ones(length_sim,1);
    time=nan*ones(length_sim, 1);
    for t=1:length_sim 
        if SKIP_AMR && mod(t*md.settings.output_frequency, md.transient.amr_frequency)==0
            continue;
        end
        sol = md.results.TransientSolution(t);
        H = sol.Thickness; %H=H(x>=CUTOFF_X);
        u = sol.Vx;
            GZ_idx = find(abs(sol.MaskOceanLevelset)<=500);
            mean_Hg(t) = mean(H(GZ_idx));
            mean_ug(t) = mean(u(GZ_idx));
        if t==1 
            mean_H(t)=mean(H);
            time(t)=md.results.TransientSolution(t).time - md.results.TransientSolution(1).time + md.timestepping.time_step;
            continue;
        end
        time(t)=md.results.TransientSolution(t).time - md.results.TransientSolution(1).time + md.timestepping.time_step;
        mean_H(t) = mean(H); 
        MS_diff_H(t) = abs(mean_H(t)- mean_H(t-1)) / (time(t) - time(t-1));
    end
    time=time(~isnan(MS_diff_H));
    MS_diff_H=MS_diff_H(~isnan(MS_diff_H));
    MS_diff_H = MS_diff_H / mean(md.smb.mass_balance);
    disp(size(MS_diff_H))

    % ---- verdict (same arithmetic the Python gate uses) --------------------
    % Two-part: terminal level at or below Threshold AND a non-increasing trend
    % over the last NUM_YEARS. Reported, never enforced -- this function does not
    % error out and does not touch any file.
    if ~isempty(THRESHOLD) && numel(MS_diff_H)>=1
        terminal = MS_diff_H(end);
        w = time >= (max(time) - NUM_YEARS);
        slope = 0;
        if sum(w) > 2
            pf = polyfit(time(w), MS_diff_H(w), 1);
            slope = pf(1);
        end
        pass = (terminal <= THRESHOLD) && (slope <= 0);
        if pass, verdict='PASS'; else, verdict='FAIL'; end
        fprintf(['[check_steady_state] :: terminal dh/dt / adot = %.4e ' ...
                 '(threshold %.1e), trend over last %g yr = %+.4e /yr => %s\n'], ...
                terminal, THRESHOLD, NUM_YEARS, slope, verdict);
    end

    if ~DO_PLOT
        return;
    end

    figure; tiledlayout(1,1);
    nexttile;
    plot(time, MS_diff_H, 'k-', 'markersize', 15, 'LineWidth', 4); 
    xlabel("Time [yr]"); title("$\partial_t h / \dot a$");
    set(gca,'YScale','log');  grid on;  % not yscale(): that needs R2023b+,
                                         % and `ml load matlab` gives R2022b

    % nexttile; hold on;
    % plot(qg_time(2:end), d_gl_flux, ...
    %             'k--', 'markersize', 15, 'LineWidth', 0.5); 
    % plot(qg_time(2:end), smoothdata(d_gl_flux,'gaussian',10), ...
    %             'k.-', 'markersize', 15, 'LineWidth', 3); 
    % title("[Mean($Q_t$) - Mean($Q_{t-1}$)]/($u_g H_g$)");
    % xlabel("Time [yr]"); yscale('log');  grid on; 

    % sgtitle("STEADY STATE CHECK", 'interpreter', 'latex', 'Color', 'red', 'FontSize', 19, 'FontWeight','bold');
    
    nexttile(1); xlim([nanmax(time)-NUM_YEARS, nanmax(time)]);
    % nexttile(2); xlim([nanmax(qg_time)-NUM_YEARS, nanmax(qg_time)]);
    set(gcf,'Position', [300 300 400 400]);
end