function md_2d_reinit = reinit_and_cut_2D(md, config)
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
    config.SAVE_PATH = [fullfile(data_root, 'cut_exp') filesep];
    config.EXECUTION_PATH = [exec_root filesep];

    disp('DEBUG: transientrestart-ing model simulation.')
    md_2d_reinit = transientrestart(md);

    md_2d_reinit.basalforcings = basalforcingsideal;
    md_2d_reinit.basalforcings.groundedice_melting_rate = zeros(md_2d_reinit.mesh.numberofvertices, 1);
    md_2d_reinit.basalforcings.geothermalflux = zeros(md_2d_reinit.mesh.numberofelements, 1);
    md_2d_reinit.basalforcings.alpha = config.alpha; 
    md_2d_reinit.basalforcings.beta  = config.beta; 
    md_2d_reinit.basalforcings.gamma = config.gamma; 
    md_2d_reinit.basalforcings.m0    = 0.0; 
    md_2d_reinit.basalforcings.width = max(md_2d_reinit.mesh.y) - min(md_2d_reinit.mesh.y);

    md_2d_reinit.transient.amr_frequency      = config.amr_config.amr_frequency;
    
    md_2d_reinit.amr.hmin                     = config.amr_config.hmin; 
    md_2d_reinit.amr.hmax                     = config.amr_config.hmax; 
    md_2d_reinit.amr.groundingline_resolution = config.amr_config.groundingline_resolution; 
    md_2d_reinit.amr.groundingline_distance   = config.amr_config.groundingline_distance;
    md_2d_reinit.amr.icefront_resolution      = config.amr_config.icefront_resolution;
    md_2d_reinit.amr.icefront_distance        = config.amr_config.icefront_distance; 

    md_2d_reinit.amr.deviatoricerror_maximum   = 500e3;
    md_2d_reinit.amr.deviatoricerror_threshold = 1.0;
    md_2d_reinit.amr.deviatoricerror_groupthreshold = 1.0;
    md_2d_reinit.amr.deviatoricerror_resolution = 1000;
    md_2d_reinit.amr.thicknesserror_maximum   = 100;
    md_2d_reinit.amr.thicknesserror_threshold = 1.0;
    md_2d_reinit.amr.thicknesserror_groupthreshold = 1.0;
    md_2d_reinit.amr.thicknesserror_resolution = 1000;
    md_2d_reinit.amr.restart = 0;
    md_2d_reinit.amr.gradation = 1.2;
    md_2d_reinit.amr.fieldname = 'None';

    icedivide_pos = find(md_2d_reinit.mesh.x <= 0.1 & md_2d_reinit.mesh.x >= -0.1);
    boundary_pos = find((md_2d_reinit.mesh.y <= min(md_2d_reinit.mesh.y)+0.5 & md_2d_reinit.mesh.y >= min(md_2d_reinit.mesh.y)-0.5) ... 
                         | (md_2d_reinit.mesh.y <= max(md_2d_reinit.mesh.y)+0.5 & md_2d_reinit.mesh.y >= max(md_2d_reinit.mesh.y)-0.5));
    md_2d_reinit.stressbalance.spcvz(icedivide_pos)=0.0;
    md_2d_reinit.stressbalance.spcvy(icedivide_pos)=0.0;
    md_2d_reinit.stressbalance.spcvx(icedivide_pos)=0.0;
    md_2d_reinit.stressbalance.spcvy(boundary_pos)=0.0;
    md_2d_reinit.stressbalance.spcvx(boundary_pos)=0.0;

    md_2d_reinit.timestepping.time_step = config.time_step;
    md_2d_reinit.timestepping.final_time = md_2d_reinit.timestepping.start_time + config.num_yrs;
    md_2d_reinit.settings.output_frequency = config.output_freq;

    md_2d_reinit.masstransport.stabilization = 5;
    md_2d_reinit.levelset.stabilization      = 5;

    md_2d_reinit.cluster = generic('name', oshostname, 'np', config.num_proc, 'executionpath', config.EXECUTION_PATH);
    req_outputs = {'default', 'StrainRatexx', 'StrainRateyy', ...
                   'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy', ...
                   'Bed', 'FrictionAlpha2'};
                   
    md_2d_reinit.transient.requested_outputs=req_outputs;
    md_2d_reinit.miscellaneous.name = [config.timestamp 'flatbed_icecut_exp' config.exp_id config.name];

    md_2d_reinit.mask.ice_levelset(md.mesh.x > config.ice_cut_x) = 1.0;
    md_2d_reinit.mask.ice_levelset = reinitializelevelset(md_2d_reinit, md_2d_reinit.mask.ice_levelset);

    md_2d_reinit=solve(md_2d_reinit,'transient');

    TOTAL_SAVE_PATH = [config.SAVE_PATH config.exp_setup '/exp_' config.exp_id '/' config.name];
    % Safeguard against very large model files
    info = whos('md_2d_reinit');
    SAVE_FILENAME = [TOTAL_SAVE_PATH '/' config.timestamp '_md_' config.name '.mat'];
    disp(['Saving results of experiment ' config.name ' to ' SAVE_FILENAME ' ...'])
    if ~exist(TOTAL_SAVE_PATH, 'dir')
        system(['mkdir -p ' TOTAL_SAVE_PATH]);
    end
    if info.bytes > 2e9
        save(SAVE_FILENAME, 'md_2d_reinit', '-v7.3');
    else
        save(SAVE_FILENAME, 'md_2d_reinit');
    end
    disp('Saved')
end