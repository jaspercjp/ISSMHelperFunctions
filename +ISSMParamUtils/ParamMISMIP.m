%PARAMMISMIP Initializes ISSM model parameters and geometry for MISMIP+ experiments.
%
%   md = PARAMMISMIP(md)
%   md = PARAMMISMIP(md, 'NoSlip', false, 'InitGeometry', true, 'StartupModel', [], 'InitPhysics', true)
%   md = PARAMMISMIP(md, 'Config', config, 'NoSlip', true)
%
%   This function sets up the ISSM model structure for MISMIP+ experiments, including
%   geometry, physical parameters, boundary conditions, and initial fields. It supports
%   initializing geometry from scratch or from a previous model, and allows configuration
%   of no-slip boundary conditions and physics initialization.
%
%   Inputs:
%       md            - ISSM model structure to be initialized.
%       NoSlip        - (Optional) Logical flag to apply no-slip at lateral walls (default: false).
%       InitGeometry  - (Optional) Logical flag to initialize geometry (default: true).
%       StartupModel  - (Optional) Previous ISSM model to initialize geometry from (default: []).
%       InitPhysics   - (Optional) Logical flag to initialize physical parameters (default: true).
%       Config        - (Optional) Experiment config struct (as produced by jsondecode of a
%                       STEADY_STATE_CONFIGS json). Fields read here:
%                         FRICTION_LAW          - 'frictionweertman' (default) | 'frictiontsai'
%                                                 | 'frictionschoof' | 'frictionregcoulomb'
%                         SURFACE_MASS_BALANCE  - [m/yr], default 0.3
%
%   Outputs:
%       md            - ISSM model structure with updated parameters and geometry.
%
%   The function also defines the bed topography using the Cornford et al. (2020) formula,
%   sets up surface mass balance, rheology, friction, and boundary conditions for the model.
%
%   Example:
%       md = ISSMParamUtils.ParamMISMIP(md, 'NoSlip', true, 'InitGeometry', false)
function md = ParamMISMIP(md, varargin)
    p = inputParser;
    addParameter(p, 'NoSlip', false, @(x) islogical(x));
    addParameter(p, 'InitGeometry', true, @(x) islogical(x));
    addParameter(p, 'InitPhysics', true, @(x) islogical(x));
    addParameter(p, 'IceTongueX', -1, @(x) isnumeric(x));
    addParameter(p, 'StartupModel', [], @(x) isempty(x) || isa(x, 'model'))
    addParameter(p, 'Config', struct(), @(x) isstruct(x));
    parse(p, varargin{:});
    NO_SLIP = p.Results.NoSlip;
    INIT_GEOMETRY=p.Results.InitGeometry; % whether
    INIT_PHYSICS=p.Results.InitPhysics;
    ICE_TONGUE_X=p.Results.IceTongueX;
    STARTUP_MODEL=p.Results.StartupModel;
    CONFIG=p.Results.Config;

    % =========== My constants =====================
    % Surface mass balance [m/yr]. Read from the experiment config so that the
    % value actually applied matches the one the config declares (and hence the
    % smb%dcmyr string calc_steady_state.m puts in the output folder name).
    % Mirrors ParamConstSlope.m. Historically SMB was hardcoded here while
    % config.SURFACE_MASS_BALANCE fed only the folder name, so the two could
    % silently disagree -- 20260615_{1737,2311}_MISMIP_frictionregcoulomb_steady_smb30cmyr
    % were both run at 0.5 m/yr despite being named smb30cmyr.
    if isfield(CONFIG, 'SURFACE_MASS_BALANCE')
        SMB=CONFIG.SURFACE_MASS_BALANCE;
    else
        SMB=0.3;    % default, matching every MISMIP config in STEADY_STATE_CONFIGS
    end
    n=3; B=1.164e8; % rheology;
    m=3; C2=3.16e6; % friction
    
    % =================== PRESCRIBE BED TOPO ====================
    md.geometry.bed = define_bed(md.mesh.x, md.mesh.y);
    di = md.materials.rho_ice / md.materials.rho_water; 

    % ================== INITIALIZE GEOMETRY ====================
    if INIT_GEOMETRY && isempty(STARTUP_MODEL)
        hmin=300;
        hmax=3000;
        xmin=min(md.mesh.x);
        xmax=max(md.mesh.x);
        md.geometry.thickness = hmax+(hmin-hmax)*(md.mesh.x-xmin)/(xmax-xmin);    
        md = ISSMModelUtils.adjust_thk_ols(md);
    elseif INIT_GEOMETRY && isa(STARTUP_MODEL, 'model')
        sol = STARTUP_MODEL.results.TransientSolution(end);
        startup_x = sol.MeshX; startup_y = sol.MeshY; 
        thk_sip = scatteredInterpolant(startup_x, startup_y, sol.Thickness, 'natural', 'nearest');
        md.geometry.thickness = thk_sip(md.mesh.x, md.mesh.y);
        md = ISSMModelUtils.adjust_thk_ols(md);
    end

    if INIT_PHYSICS
        % ---------------- friction parametrization (config-selected) ----------------
        if isfield(CONFIG, 'FRICTION_LAW')
            friction_law = CONFIG.FRICTION_LAW;
        else
            disp('[ParamMISMIP] :: No friction parametrization specified. Defaulting to Weertman.')
            friction_law = 'frictionweertman';
        end

        if strcmp(friction_law, 'frictionweertman')
            md.friction = frictionweertman();
            m = 3.0;
            coeff = sqrt(C2);
        elseif strcmp(friction_law, 'frictiontsai')
            md.friction = frictiontsai();
            md.mask.ocean_levelset = reinitializelevelset(md, md.mask.ocean_levelset);
            md.friction.f = 0.8 * ones(md.mesh.numberofvertices, 1);
            md.friction.u0 = 600;
            md.friction.groundingline_coulomb_distance = 15e3;
            m = 3.0; % m is defined differently for frictiontsai
            coeff = sqrt(C2); % C is not squared in the tsai parametrization
        elseif strcmp(friction_law, 'frictionschoof')
            md.friction = frictionschoof();
            md.friction.Cmax = 0.6 * ones(md.mesh.numberofvertices, 1);
            m = 1/3.0;
            coeff = sqrt(C2);
        elseif strcmp(friction_law, 'frictionregcoulomb')
            md.friction = frictionregcoulomb();
            m = 3.0;
            coeff = sqrt(C2);
            md.friction.u0 = 200;
        else
            error(['Unsupported friction parametrization ' friction_law])
        end
        % ----------------------------------------------------------------------------
        md.materials.rheology_n = nan*ones(md.mesh.numberofelements, 1);
        md.materials.rheology_B = nan*ones(md.mesh.numberofvertices, 1);
        md.friction.m =           nan*ones(md.mesh.numberofelements, 1);
        md.friction.C =           nan*ones(md.mesh.numberofvertices, 1);

        md.friction.m(:)           = m;
        md.friction.C(:)           = coeff;
        md.materials.rheology_B(:) = B;
        md.materials.rheology_n(:) = n;
        
        md=ISSMParamUtils.init_BC(md, 'NoSlip', NO_SLIP); % Note this is free-slip by default

        % boundary conditions for mass balance 
        md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
        md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
        md.groundingline.migration='SubelementMigration';
        md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
        md.smb.mass_balance=SMB*ones(md.mesh.numberofvertices, 1);
    end

    md.transient.issmb = 1;
    md.transient.ismasstransport = 1;
    md.transient.isthermal = 0;
    md.transient.isstressbalance = 1;
    md.transient.isgroundingline = 1;
    md.transient.ismovingfront = 0;
end

function [bed] = define_bed(x,y)
    % ========= Constants from Cornford et al. (2020)
    B0=-150; B2=-728.8; B4=343.91; B6=-50.75;
    wc=24000; fc=4000; dc=500;
    xbar=300000;
    % ==========
    md_x = x/xbar;
    Bx = B0 + B2*md_x.^2 + B4*md_x.^4 + B6*md_x.^6;
    By = dc*(1./(1+exp(-2*(y-wc)/fc)) + 1./(1+exp(2*(y+wc)/fc)));
    bed=max(Bx+By, -720);
end