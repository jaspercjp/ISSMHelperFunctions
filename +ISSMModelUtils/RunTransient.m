%% Use this function to setup and run a transient simulation of an ice shelf
function md=RunTransient(model_name, domain_file, param_file, varargin)
    % Handle named inputs 
    model_name=char(model_name);
    domain_file=char(domain_file);
    param_file=char(param_file);
    p = inputParser;

    % model setup
    addParameter(p, 'NumProc', 1, @(x) isnumeric(x) && x>0);
    addParameter(p, 'SaveModel', true, @(x) islogical(x));
    addParameter(p, 'RestartModel', [], @(x) isempty(x) || isa(x,'model'));
    addParameter(p, 'SavePrevious', false, @(x) islogical(x));
    addParameter(p, 'ExecutionPath', '', @(x) true);
    addParameter(p, 'InterpolateForcings', false, @(x) islogical(x));

    % specific model settings 
    addParameter(p, 'UseAMR', false, @(x) islogical(x));
    addParameter(p, 'AMRFrequency', 0, @(x) isnumeric(x));
    addParameter(p, 'AMRSetting', [], @(x) isempty(x) || isstruct(x))
    addParameter(p, 'Hmax', 10000, @(x) isnumeric(x) && x>0);
    addParameter(p, 'Hmin', 1000, @(x) isnumeric(x) && x>0);

    addParameter(p, 'MeltRate', 0.0, @(x) isnumeric(x) || isnan(x));
    addParameter(p, 'DeepMeltRate', nan, @(x) isnumeric(x)||isnan(x));
    addParameter(p, 'QuadMelt', [], @(x) isa(x, 'basalforcingsismip6')||isempty(x));
    addParameter(p, 'SurfaceMassBalance', nan, @(x) isnumeric(x)||isnan(x));

    addParameter(p, 'TimeLength', 200, @(x) isnumeric(x) && x>0);
    addParameter(p, 'TimeStep', 0.2, @(x) isnumeric(x) && x>0);
    addParameter(p, 'OutputFreq', 50, @(x) isnumeric(x) && x>0);
    addParameter(p, 'SavePath', "../saved/", @(x) isstring(x) || ischar(x));

    % physics components
    addParameter(p, 'Calving', "NONE", @(x) isstring(x) || ischar(x));
    addParameter(p, 'CalvingConfig', [], @(x) isempty(x) || isstruct(x));
    addParameter(p, 'IsThermal', false, @(x) islogical(x));
    addParameter(p, 'ThermalConfig', [], @(x) isempty(x) || isstruct(x));
    parse(p, varargin{:});
    
    % Extract Parameters
    NUM_PROC = p.Results.NumProc;
    SAVE_MODEL = p.Results.SaveModel;
    HMAX = p.Results.Hmax;
    MELT_RATE = p.Results.MeltRate;
    DEEP_MELT_RATE = p.Results.DeepMeltRate;
    QUAD_MELT = p.Results.QuadMelt;
    TIME_LENGTH = p.Results.TimeLength;
    OUTPUT_FREQ = p.Results.OutputFreq;
    TIME_STEP = p.Results.TimeStep;
    SAVE_PATH = p.Results.SavePath;
    CALVING = p.Results.Calving;
    HMIN = p.Results.Hmin;
    USE_AMR = p.Results.UseAMR;
    AMR_FREQ=p.Results.AMRFrequency;
    AMR_SETTING=p.Results.AMRSetting;
    RESTART_MODEL=p.Results.RestartModel;
    CALVING_CONFIG=p.Results.CalvingConfig;
    SURFACE_MASS_BALANCE=p.Results.SurfaceMassBalance;
    INTERPOLATE_FORCINGS=p.Results.InterpolateForcings;
    SAVE_PREVIOUS = p.Results.SavePrevious;
    EXECUTION_PATH = p.Results.ExecutionPath;
    IS_THERMAL = p.Results.IsThermal;
    THERMAL_CONFIG = p.Results.ThermalConfig;

    % Setup. Handle the case where we want to restart the run from an
    % existing model
    if ~isa(RESTART_MODEL,'model')
        md=bamg(model,'domain',domain_file,'hmax',HMAX,'hmin',HMIN,'splitcorners',1);
        md=parameterize(md,param_file);
        md=setflowequation(md,'SSA','all');
        % Transient settings
        md.settings.output_frequency=OUTPUT_FREQ; 
        md.timestepping.final_time = TIME_LENGTH;
        md.timestepping.time_step = TIME_STEP;
    else
        md=RESTART_MODEL;
        md.settings.output_frequency=OUTPUT_FREQ; 
        md.timestepping.final_time = md.timestepping.final_time+TIME_LENGTH;
        md.timestepping.time_step = TIME_STEP;
    end

    % Giving each job an unique task name is important if we are running 
    % multiple jobs at the same time
    rng('shuffle');
    split_str_tmp = split(SAVE_PATH,'/');
    md.miscellaneous.name = [char(split_str_tmp{end-1})  '_rndi_' num2str(randi(50000))];
    
    % Stressbalance settings
    md.stressbalance.maxiter=30;
	% md.stressbalance.abstol=NaN;
	% md.stressbalance.restol=1;
    % md.settings.solver_residue_threshold=1e-2; % we are relaxing the solver residue threshold here

    % Automatic mesh refinement (AMR) settings
    if USE_AMR && isstruct(AMR_SETTING)
        md.transient.amr_frequency=AMR_SETTING.amr_freq;
        disp('[RunTransient] :: AMR Settings provided. Reading in options...')
        setting_name=fieldnames(AMR_SETTING);
        for i=1:length(setting_name)
            if ~strcmp(setting_name{i}, 'amr_freq')
                md.amr.(setting_name{i}) = AMR_SETTING.(setting_name{i});
            end
        end
    elseif USE_AMR && ~isa(AMR_SETTING, 'amr')
        % use default AMR settings 
        disp('[RunTransient] :: No AMR options specified -- Using default AMR settings.')
        md.transient.amr_frequency=1;
        md.amr.groundingline_resolution=2000;
        md.amr.groundingline_distance=20000;
        md.amr.icefront_resolution=100;
        md.amr.icefront_distance=5000;
        md.amr.hmin=HMIN; % the same resolution used around the grounding line
        md.amr.hmax=HMAX; % the same coase resolution used to generate the coarse mesh
        md.amr.gradation=3; % this controls the ratio between two consecutive edges
        md.amr.fieldname='None'; 
        md.amr.keepmetric=0; % no field, no metric
    elseif ~USE_AMR
        md.transient.amr_frequency=0;
    end

    % Important to set the "name" of the simulation for ISSM to register its 
    % executable path
    if isempty(EXECUTION_PATH)
        md.cluster = generic('name', oshostname(), 'np', NUM_PROC);
    else
        md.cluster=generic('name',oshostname(),'np', NUM_PROC, 'executionpath', EXECUTION_PATH);
    end
    md.transient.requested_outputs = {'default', 'Bed', 'MaskIceLevelset',...
                                            'DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy',...
                                            'StrainRatexx','StrainRatexy','StrainRateyy',...
                                            'FloatingArea', ...
                                            'BasalforcingsFloatingiceMeltingRate'};
    
    
    if INTERPOLATE_FORCINGS
        md.timestepping.interp_forcing = 1;
    else
        md.timestepping.interp_forcing = 0;
    end
    
    % Specify basal melt 
    if isscalar(MELT_RATE) && isnan(DEEP_MELT_RATE) && ~isnan(MELT_RATE)
        % use spatially uniform melt 
        md.basalforcings.floatingice_melting_rate=MELT_RATE*ones(md.mesh.numberofvertices, 1);
    elseif isscalar(MELT_RATE) && ~isnan(DEEP_MELT_RATE)
        % use linear-depth-dependence melt 
        md.basalforcings = linearbasalforcings();
        md.basalforcings.deepwater_melting_rate = DEEP_MELT_RATE;
        md.basalforcings.upperwater_melting_rate = MELT_RATE;
        md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);
        md.deepwater_elevation = -250;
        md.upperwater_elevation = -100;
    elseif all(isnan(MELT_RATE)) && isnan(DEEP_MELT_RATE) && all(class(QUAD_MELT)=='basalforcingsismip6')
        disp("- [RunTransient] :: Using quadratic melt parametrization.")
        % non-local quadratic parametrization
        md.basalforcings = QUAD_MELT;
    elseif length(MELT_RATE)==md.mesh.numberofvertices 
        % Spatially varying melt rate 
        md.basalforcings.floatingice_melting_rate=MELT_RATE;
    else
        % Time dependent melt rate... 
        % NOTE: This approach is currently NOT supported by ISSM 
        % In this case, the user is assumed to have passed in a forcing array with
        % the shape (2, # of forcing values), the first row being the
        % spatially uniform melt rate, and the second being times at which
        % these forcings are applied
        md.basalforcings.floatingice_melting_rate = repmat(MELT_RATE(1,:), [md.mesh.numberofvertices, 1]);
        md.basalforcings.floatingice_melting_rate(end+1,:) = MELT_RATE(2,:);
        md.transient.requested_outputs = [md.transient.requested_outputs, {'BasalforcingsFloatingiceMeltingRate'}];
    end

    if IS_THERMAL
        md.transient.isthermal = 1;
        md.thermal.spctemperature = nan*ones(md.mesh.numberofvertices, 1);
        md.basalforcings.geothermalflux = THERMAL_CONFIG.geothermalflux * ones(md.mesh.numberofvertices, 1);
        md.initialization.temperature = THERMAL_CONFIG.initial_temperature * ones(md.mesh.numberofvertices, 1);
        md.initialization.waterfraction = zeros(md.mesh.numberofvertices, 1);
    end 

    if ~isnan(SURFACE_MASS_BALANCE)
        md.smb.mass_balance = SURFACE_MASS_BALANCE * ones(md.mesh.numberofvertices, 1);
    end
    
    % Specify calving parametrization
    if ~(CALVING=="NONE" || isempty(CALVING))
        md.transient.ismovingfront = 1;
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        if CALVING=="VM"
            md.calving=calvingvonmises();
            if ~isempty(CALVING_CONFIG)
                md.calving.stress_threshold_floatingice = CALVING_CONFIG.stress_threshold_floatingice;
                md.calving.min_thickness = CALVING_CONFIG.min_thickness;
            else
                md.calving.stress_threshold_floatingice = 150e3; % 150 kPa is the value for Amery 
                md.calving.min_thickness = 50;
            end
            md.levelset.stabilization = 1;
            md.levelset.reinit_frequency=10;
            md.levelset.migration_max = 10e9;
        elseif CALVING=="EC"
            md.calving=calvinglevermann();
            if ~isempty(CALVING_CONFIG)
                md.calving.coeff = CALVING_CONFIG.coeff*md.constants.yts *ones(md.mesh.numberofvertices,1); 
            else
                md.calving.coeff = 5e7*md.constants.yts *ones(md.mesh.numberofvertices,1); 
            end
            md.levelset.stabilization = 1;
            md.levelset.reinit_frequency=10;
            md.levelset.migration_max = 10e9;
        elseif CALVING=="BUTTRESS"
            md.calving=calvingcrevassedepth();
            if ~isempty(CALVING_CONFIG)
                md.calving.crevasse_threshold = CALVING_CONFIG.crevasse_threshold;
                md.calving.crevasse_opening_stress = CALVING_CONFIG.crevasse_opening_stress;
                md.calving.water_height=CALVING_CONFIG.water_height*ones(md.mesh.numberofvertices,1); 
                md.calving.timescale=CALVING_CONFIG.timescale;
                if isfield(CALVING_CONFIG, 'min_iceberg_size')
                    md.calving.min_iceberg_size=CALVING_CONFIG.min_iceberg_size;
                else
                    md.calving.min_iceberg_size=5e3;
                end
                if isfield(CALVING_CONFIG, 'advect_icefront')
                    md.calving.advect_icefront=CALVING_CONFIG.advect_icefront;
                else
                    md.calving.advect_icefront=0;
                end
            else
                md.calving.crevasse_threshold = 1.0; 
                md.calving.crevasse_opening_stress=3;
                md.calving.timescale=-1;
                % [TODO]: Incorporate water depth into the ISSM code in Tria.cpp 
                md.calving.water_height=zeros(md.mesh.numberofvertices,1); 
            end
            % [TODO]: What are the best levelset simulation parameters to use?
            md.levelset.stabilization = 1;
            md.levelset.reinit_frequency=floor(md.calving.timescale / md.timestepping.time_step);
            md.transient.requested_outputs=[md.transient.requested_outputs, {'CrevasseDepth', 'ButtressingK', 'CalvingPropagatedArea', 'CalvingPropagatedMinX'}];
        elseif CALVING=="MINTHICKNESS"
            md.calving=calvingminthickness();
            if ~isempty(CALVING_CONFIG)
                md.calving.min_thickness=CALVING_CONFIG.min_thickness;
            else
                md.calving.min_thickness=80;
            end
        elseif CALVING=="MOVEFRONT"
            % Just advect the ice front without calving.
            md.levelset.reinit_frequency=10;
            md.calving.calvingrate = zeros(md.mesh.numberofvertices, 1);
        else
            error("Unsupported Calving Parametrization.");
        end
    else
        disp('- [RunTransient] :: Empty calving specification, defaulting to fixed calving front.')
        md.transient.ismovingfront = 0;
    end

    % Solve the transient problem
    disp('- [RunTransient] :: Solving transient...')
    % md.verbose.solution = false;  
    % md.verbose.module = false;  
    % md.verbose.mprocessor = false;  
    md=solve(md,'Transient');
    

    if length(fieldnames(md.results))>1
        if isa(RESTART_MODEL, 'model') && ~SAVE_PREVIOUS
            % Don't save redundant data if we are only restarting the run
            allFields = fieldnames(md.results);
            fieldToKeep = 'TransientSolution';
            fieldsToDelete = allFields(~strcmp(allFields, fieldToKeep));
            md.results = rmfield(md.results, fieldsToDelete);
        elseif SAVE_PREVIOUS
            % 1. Get all TransientSolution field names
            allFields = fieldnames(md.results);
            pattern = 'TransientSolution\d*$';
            transientFields = allFields(~cellfun('isempty', regexp(allFields, pattern, 'once')));

            % 2. Get indices and sort (0 for base, N for numbered)
            fieldIndices = zeros(1, numel(transientFields));
            for i = 1:numel(transientFields)
                tokens = regexp(transientFields{i}, 'TransientSolution(\d*)', 'tokens', 'once');
                if ~isempty(tokens{1})
                    fieldIndices(i) = str2double(tokens{1});
                end
            end
            [~, sortIdx] = sort(fieldIndices);

            % 3. Reorder: Move the base 'TransientSolution' (index 0) to the end
            sortedFields = transientFields(sortIdx);
            idxBase = find(fieldIndices(sortIdx) == 0, 1);
            reorderedFields = [sortedFields(idxBase+1:end); sortedFields(1:idxBase)];

            % 4. Concatenate
            finalStruct = {};
            for name = reorderedFields'
                finalStruct = [finalStruct, {md.results.(name{1})}];
            end

            md.results = rmfield(md.results, transientFields); % Clean up original fields
            md.results = rmfield(md.results, 'reinitialization');
            md.results.TransientSolution = finalStruct; % Result stored here
        end
    end
    
    if SAVE_MODEL
        if ~exist(SAVE_PATH, 'dir')
            system(['mkdir -p ' SAVE_PATH]);
        end
        disp(['- [RunTransient] :: Saving to ' SAVE_PATH model_name]);
        info = whos('md');
        % Safeguard against very large model files
        if info.bytes > 2e9
            save([SAVE_PATH model_name], 'md', '-v7.3');
        else
            save([SAVE_PATH model_name], 'md');
        end
    end
end