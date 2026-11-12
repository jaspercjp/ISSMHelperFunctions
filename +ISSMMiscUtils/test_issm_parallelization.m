function results = test_issm_parallelization(md, processor_counts, num_trials, varargin)
% TEST_ISSM_PARALLELIZATION - Test ISSM stress balance performance with different processor counts
%
% USAGE:
%   results = test_issm_parallelization(md, processor_counts, num_trials)
%   results = test_issm_parallelization(md, processor_counts, num_trials, 'cluster', cluster_obj)
%
% INPUTS:
%   md              - ISSM model object (already parametrized)
%   processor_counts - Vector of processor counts to test [1, 2, 4, 8, ...]
%   num_trials      - Number of trials per processor count for statistics
%   
% OPTIONAL:
%   'cluster'       - Custom cluster object (default: generic cluster)
%   'verbose'       - Display progress (default: true)
%   'save_results'  - Save results to file (default: false)
%   'output_file'   - Output filename (default: 'issm_parallel_results.mat')
%
% OUTPUTS:
%   results - Structure with timing statistics and analysis
%
% EXAMPLE:
%   proc_counts = [1, 2, 4, 8, 16];
%   results = test_issm_parallelization(md, proc_counts, 5);
%   plot_parallel_results(results);

% Parse input arguments
p = inputParser;
addRequired(p, 'md');
addRequired(p, 'processor_counts', @(x) isnumeric(x) && all(x > 0));
addRequired(p, 'num_trials', @(x) isnumeric(x) && x > 0);
addParameter(p, 'cluster', [], @(x) isa(x, 'cluster'));
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'save_results', false, @islogical);
addParameter(p, 'output_file', 'issm_parallel_results.mat', @ischar);
parse(p, md, processor_counts, num_trials, varargin{:});

% Extract parsed values
use_custom_cluster = ~isempty(p.Results.cluster);
verbose = p.Results.verbose;
save_results = p.Results.save_results;
output_file = p.Results.output_file;

% Initialize results structure
n_proc_tests = length(processor_counts);
results = struct();
results.processor_counts = processor_counts;
results.num_trials = num_trials;
results.wall_times = cell(n_proc_tests, 1);
results.mean_times = zeros(n_proc_tests, 1);
results.std_times = zeros(n_proc_tests, 1);
results.min_times = zeros(n_proc_tests, 1);
results.max_times = zeros(n_proc_tests, 1);
results.speedup = zeros(n_proc_tests, 1);
results.efficiency = zeros(n_proc_tests, 1);
results.test_date = datestr(now);

if verbose
    fprintf('\n=== ISSM Stress Balance Parallelization Test ===\n');
    fprintf('Model elements: %d\n', md.mesh.numberofelements);
    fprintf('Model vertices: %d\n', md.mesh.numberofvertices);
    fprintf('Processor counts to test: [%s]\n', num2str(processor_counts));
    fprintf('Trials per configuration: %d\n', num_trials);
    fprintf('================================================\n\n');
end

% Store original cluster settings
original_cluster = md.cluster;

% Test each processor count
for i = 1:n_proc_tests
    nproc = processor_counts(i);
    trial_times = zeros(num_trials, 1);
    
    if verbose
        fprintf('Testing %d processors (%d/%d)...\n', nproc, i, n_proc_tests);
    end
    
    % Configure cluster for current processor count
    if use_custom_cluster
        md.cluster = p.Results.cluster;
        md.cluster.np = nproc;
    else
        md.cluster = generic('name', oshostname(), 'np', nproc);
    end
    
    % Run multiple trials
    for trial = 1:num_trials
        if verbose && num_trials > 1
            fprintf('  Trial %d/%d... ', trial, num_trials);
        end
        
        % Create a copy of the model to avoid any state changes
        md_trial = md;
        % md_trial.verbose=verbose(false);
        % md_trial.verbose.solution=true;
        
        % Time the stress balance solve
        tic;
        md_trial = solve(md_trial, 'stressbalance');
        trial_time = toc;
        
        trial_times(trial) = trial_time;
        
        if verbose && num_trials > 1
            fprintf('%.2f s\n', trial_time);
        end
    end
    
    % Calculate statistics
    results.wall_times{i} = trial_times;
    results.mean_times(i) = mean(trial_times);
    results.std_times(i) = std(trial_times);
    results.min_times(i) = min(trial_times);
    results.max_times(i) = max(trial_times);
    
    if verbose
        fprintf('  Mean: %.2f ± %.2f s (min: %.2f, max: %.2f)\n', ...
               results.mean_times(i), results.std_times(i), ...
               results.min_times(i), results.max_times(i));
    end
end

% Calculate speedup and efficiency relative to serial case
serial_time = results.mean_times(1);  % Assume first entry is serial (1 proc)
for i = 1:n_proc_tests
    results.speedup(i) = serial_time / results.mean_times(i);
    results.efficiency(i) = results.speedup(i) / processor_counts(i) * 100;
end

% Restore original cluster
md.cluster = original_cluster;

% Display summary
if verbose
    fprintf('\n=== PERFORMANCE SUMMARY ===\n');
    fprintf('Procs | Mean Time (s) | Std Dev (s) | Speedup | Efficiency (%%)\n');
    fprintf('------|---------------|-------------|---------|---------------\n');
    for i = 1:n_proc_tests
        fprintf('%5d | %11.2f | %9.2f | %6.2fx | %11.1f\n', ...
               processor_counts(i), results.mean_times(i), results.std_times(i), ...
               results.speedup(i), results.efficiency(i));
    end
    fprintf('===========================\n\n');
end

% Save results if requested
if save_results
    save(output_file, 'results');
    if verbose
        fprintf('Results saved to: %s\n', output_file);
    end
end

% Create summary plot
create_performance_plots(results);

end

function create_performance_plots(results)
% Create performance visualization plots

figure('Position', [100, 100, 1200, 400]);

% Plot 1: Execution time with error bars
subplot(1, 3, 1);
errorbar(results.processor_counts, results.mean_times, results.std_times, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Number of Processors');
ylabel('Execution Time (seconds)');
title('Execution Time vs Processors');
grid on;

% Plot 2: Speedup
subplot(1, 3, 2);
plot(results.processor_counts, results.speedup, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(results.processor_counts, results.processor_counts, '--k', 'LineWidth', 1); % Ideal speedup
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Number of Processors');
ylabel('Speedup');
title('Speedup vs Processors');
legend('Actual', 'Ideal', 'Location', 'northwest');
grid on;

% Plot 3: Parallel efficiency
subplot(1, 3, 3);
plot(results.processor_counts, results.efficiency, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot([min(results.processor_counts), max(results.processor_counts)], [100, 100], '--k', 'LineWidth', 1);
set(gca, 'XScale', 'log');
xlabel('Number of Processors');
ylabel('Parallel Efficiency (%)');
title('Parallel Efficiency vs Processors');
ylim([0, 110]);
grid on;

sgtitle('ISSM Stress Balance Parallelization Performance', 'FontSize', 14, 'FontWeight', 'bold');
end