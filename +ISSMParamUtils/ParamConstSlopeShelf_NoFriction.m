% slope - in radians. 
% Positive slope correspond to retrograde and negative prograde. The bed will 
% be entirely below sea level. 
function md=ParamConstSlopeShelf_NoFriction(md, hmin, hmax, slope, offset)
    disp('      defining geometry');
    
    xmin=min(md.mesh.x);
    xmax=max(md.mesh.x);
    md.geometry.thickness = hmax+(hmin-hmax)*(md.mesh.x-xmax)/(xmax-xmin);
    
    if slope < 0 
        md.geometry.bed = -tan(abs(slope)) * md.mesh.x - offset;
    else
        md.geometry.bed = tan(slope)*md.mesh.x - tan(slope)*max(md.mesh.x) - offset;
    end
    
    md.geometry.base = md.geometry.bed;
    floating_base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
    pos = find(floating_base>md.geometry.bed);
    md.geometry.base(pos) = floating_base(pos);
    md.geometry.surface = md.geometry.base + md.geometry.thickness;
    
    % Adjust the ocean level set mask
    di=md.materials.rho_ice/md.materials.rho_water;
    md.geometry.thickness=md.geometry.surface-md.geometry.base;
    md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;
    
    disp('   setting pressure');
    md.initialization.pressure=md.materials.rho_ice*md.constants.g*md.geometry.thickness;

    disp('      initial velocity');
    md.initialization.vx=zeros(md.mesh.numberofvertices,1);
    md.initialization.vy=zeros(md.mesh.numberofvertices,1);
    md.initialization.vz=zeros(md.mesh.numberofvertices,1);
    md.initialization.vel=zeros(md.mesh.numberofvertices,1);
    
    disp('      creating flow law parameter');
    md.materials.rheology_B=1e8*ones(md.mesh.numberofvertices,1); % Taken from Gudmundsson et al. (2012)
    md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
    md.damage.D=zeros(md.mesh.numberofvertices,1);
    
    disp('      creating boundary conditions');
    % Make sure the open ocean front stays an open ocean front
    md.levelset.spclevelset=nan*ones(md.mesh.numberofvertices, 1);
    md.mask.ice_levelset = -1 * ones(md.mesh.numberofvertices, 1);
    open_ocean = md.mesh.vertexonboundary & (md.mesh.x==max(md.mesh.x));
    md.mask.ice_levelset(open_ocean) = 1.0;
    md.levelset.spclevelset(open_ocean) = 1.0;
    
    % Initialize some stressbalance boundary conditions
    md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
    md.stressbalance.loadingforce=zeros(md.mesh.numberofvertices,3);
    md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
    md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
    md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
    
    % Mass transport constraints
    md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
    
    % Dirichlet on the sides up to the open ocean front
    md.stressbalance.spcvx(md.mesh.vertexonboundary & ~open_ocean) = 0;
    md.stressbalance.spcvy(md.mesh.vertexonboundary & ~open_ocean) = 0;
    md.stressbalance.spcvz(md.mesh.vertexonboundary & ~open_ocean) = 0;
    
    % Use zero basal forcing for now
    md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
    md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
    md.groundingline.migration='SubelementMigration';
    
    disp('      setting up transient simulation')
    md.transient.issmb = 1;
    md.transient.ismasstransport = 1;
    md.transient.isthermal = 0;
    md.transient.isstressbalance = 1;
    md.transient.isgroundingline = 1;
    
    md.smb.mass_balance = zeros(md.mesh.numberofvertices, 1);
end