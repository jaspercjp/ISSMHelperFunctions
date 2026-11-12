% slope - in radians. 
% Positive slope correspond to retrograde and negative prograde. The bed will 
% be entirely below sea level. 
function md=ParamConstSlope(md, slope, varargin)
    p = inputParser;
    % addParameter(p, 'NoSlip', false, @(x) islogical(x));
    addParameter(p, 'InitGeometry', true, @(x) islogical(x));
    addParameter(p, 'InitPhysics', true, @(x) islogical(x));
    addParameter(p, 'Offset', 500, @(x) isnumeric(x));
    addParameter(p, 'Hmin', 300, @(x) isnumeric(x));
    addParameter(p, 'Hmax', 3000, @(x) isnumeric(x));
    addParameter(p, 'Config', struct(), @(x) isstruct(x));
    addParameter(p, 'InitSteady1DFile', '', @(x) ischar(x));
    % addParameter(p, 'IceTongueX', -1, @(x) isnumeric(x));
    parse(p, varargin{:});
    INIT_GEOMETRY = p.Results.InitGeometry;
    INIT_PHYSICS = p.Results.InitPhysics; 
    offset=p.Results.Offset;
    CONFIG=p.Results.Config;
    INIT_STEADY_1D_FILE=p.Results.InitSteady1DFile;

    hmin=p.Results.Hmin;
    hmax=p.Results.Hmax;

    % Defaults, overridable per-field from CONFIG below.
    n=3.0; B=1e8;       % rheology
    m=3.0; C2=7.6247e6; % friction

    % ================== INITIALIZE GEOMETRY ====================
    md.geometry.bed = -offset + slope * md.mesh.x ;
    di=md.materials.rho_ice/md.materials.rho_water;
    
    if INIT_GEOMETRY && isempty(INIT_STEADY_1D_FILE)
        % Initialize with a linear thickness profile if no 1D file is provided
        xmin=min(md.mesh.x);
        xmax=max(md.mesh.x);
        md.geometry.thickness = hmax+(hmin-hmax)*(md.mesh.x-xmin)/(xmax-xmin);    
        md.geometry.base = md.geometry.bed;
        floating_base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
        pos = find(floating_base>md.geometry.bed);
        md.geometry.base(pos) = floating_base(pos);
        md.geometry.surface = md.geometry.base + md.geometry.thickness; 
        md.geometry.thickness=md.geometry.surface-md.geometry.base;
        md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;
    end

    % ======================= INITIALIZE PHYSICS ================
    if INIT_PHYSICS
        % !!! If we change the rheological exponent, all material properties must be changed !!!
        % Note that the parametrization functions use n=3 as default! 
        % ---------------- friction parametrization (config-selected) ----------------
        if isfield(CONFIG, 'FRICTION_LAW')
            friction_law = CONFIG.FRICTION_LAW;
        else
            disp('[ParamConstSlope] :: No friction parametrization specified. Defaulting to Weertman.')
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

        % Per-field config overrides (independent of one another, so that a
        % config may set FRICTION_LAW without also having to set RHEOLOGY_*).
        if isfield(CONFIG, 'RHEOLOGY_N')
            n = CONFIG.RHEOLOGY_N;
            if ~strcmp(friction_law, 'frictionschoof')
                m = CONFIG.RHEOLOGY_N;
            end
        end
        if isfield(CONFIG, 'RHEOLOGY_B')
            B = CONFIG.RHEOLOGY_B;
        end
        if isfield(CONFIG, 'FRICTION_C')
            coeff = sqrt(CONFIG.FRICTION_C);
        end
        % ----------------------------------------------------------------------------

        md.materials.rheology_n = nan*ones(md.mesh.numberofelements, 1);
        md.friction.m =           nan*ones(md.mesh.numberofelements, 1);
        md.materials.rheology_B = nan*ones(md.mesh.numberofvertices, 1);
        md.friction.C =           nan*ones(md.mesh.numberofvertices, 1);
        md.friction.m(:)           = m;
        md.friction.C(:)           = coeff;
        md.materials.rheology_B(:) = B;
        md.materials.rheology_n(:) = n;

        md=ISSMParamUtils.init_BC(md, 'NoSlip', true); % No-slip by default
        
        % Mass transport constraints
        md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
        md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
        md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
        md.groundingline.migration='SubelementMigration';
        if isfield(CONFIG, 'SURFACE_MASS_BALANCE')
            md.smb.mass_balance = CONFIG.SURFACE_MASS_BALANCE * ones(md.mesh.numberofvertices, 1);
        else
            % Vector, not a scalar: the if-branch above sets a per-vertex field and
            % ISSM's checkconsistency expects [numberofvertices,1] here too.
            md.smb.mass_balance = 0.3 * ones(md.mesh.numberofvertices, 1);
        end
        
        % set up transient solution
        md.transient.issmb = 1;
        md.transient.ismasstransport = 1;
        md.transient.isthermal = 0;
        md.transient.isstressbalance = 1;
        md.transient.isgroundingline = 1;
        md.transient.ismovingfront = 0;
    end
end