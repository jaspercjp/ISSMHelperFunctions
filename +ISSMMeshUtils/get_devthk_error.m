function [eps_tau, eps_thk] = get_devthk_error(md, config)
    % Output + execution roots come from the environment (exported by the
    % launcher: $DATA_ROOT / $ISSM_EXECUTION_ROOT), and an unset
    % DATA_ROOT is an error rather than a hard-coded site path.
    data_root = getenv('DATA_ROOT');
    if isempty(data_root)
        error('ISSMHelperFunctions:noDataRoot', 'DATA_ROOT is not set; the launchers (Core.issm_run, _calc_steady_state.sh) export it.');
    end
    exec_root = getenv('ISSM_EXECUTION_ROOT');
    if isempty(exec_root)
        exec_root = fullfile(data_root, 'ISSM_execution');
    end
    config.SAVE_PATH = [fullfile(data_root, 'melt_exp') filesep];
    config.EXECUTION_PATH = [exec_root filesep];

    md = transientrestart(md); 
    md.transient.amr_frequency = 1;
    md.amr.hmin                     = config.amr_config.hmin; 
    md.amr.hmax                     = config.amr_config.hmax; 
    md.amr.groundingline_resolution = config.amr_config.groundingline_resolution; 
    md.amr.groundingline_distance   = config.amr_config.groundingline_distance;
    md.amr.icefront_resolution      = config.amr_config.icefront_resolution;
    md.amr.icefront_distance        = config.amr_config.icefront_distance; 
    md.timestepping.final_time = md.timestepping.final_time + 2*md.timestepping.time_step;
    md.cluster = generic('name', oshostname, 'np', config.num_proc, 'executionpath', config.EXECUTION_PATH);
    % The run name is the ISSM execution-directory key AND the marshalling
    % .bin/.outbin basename, so it must be unique across every job that can run
    % concurrently against the same $ISSM_EXECUTION_ROOT.
    %
    % It was ['get_devthk_error_' timestamp '_exp' exp_id name], which omits
    % exp_setup -- so flatbed's exp_2d_1c and MISMIP's exp_2d_1c collided
    % whenever they shared a timestamp, which every batch launched in the same
    % second does. Measured 2026-07-29: 9 of 30 melt jobs died in here, split
    % between "File get_devthk_error_....bin not found" (a sibling deleted it)
    % and "Unrecognized field name DeviatoricStressErrorEstimator" (read back a
    % sibling's results). exp_setup makes it unique per family; the pid makes it
    % unique even for two runs of the SAME experiment in one allocation.
    md.miscellaneous.name = ['get_devthk_error_' config.timestamp '_' config.exp_setup ...
                             '_exp' config.exp_id config.name '_p' num2str(feature('getpid'))];
    md = solve(md,'transient');
    eps_tau = max(md.results.TransientSolution(1).DeviatoricStressErrorEstimator);
    eps_thk = max(md.results.TransientSolution(1).ThicknessErrorEstimator);
end