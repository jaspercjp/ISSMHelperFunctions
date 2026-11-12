%INIT_BC Initializes boundary conditions and fields for an ISSM model.
%
%   md = INIT_BC(md, no_slip)
%
%   This function sets up initial fields and boundary conditions for an ISSM model,
%   including ice levelset, velocity, pressure, and stress balance constraints.
%   It applies no-slip or free-slip conditions at lateral walls, zero velocity at
%   the ice divide, and water pressure at the ice front.
%
%   Inputs:
%       md       - ISSM model structure to be initialized.
%       no_slip  - Logical flag to apply no-slip at lateral walls (true/false).
%
%   Outputs:
%       md       - ISSM model structure with updated boundary conditions and fields.
%
%   The function also reinitializes the ice levelset and sets special levelset values
%   at the ice
function md = init_BC(md, varargin)
    p = inputParser;
    addParameter(p, 'NoSlip', false, @(x) islogical(x));
    addParameter(p, 'IceTongueX', -1, @(x) isnumeric(x));
    parse(p, varargin{:});
    NO_SLIP = p.Results.NoSlip;
    ICE_TONGUE_X=p.Results.IceTongueX;

    % warning('[init_BC] :: Assumes that the longitudinal direction is the x-axis!')
    md.mask.ice_levelset=-1*ones(md.mesh.numberofvertices, 1);
    md.levelset.spclevelset=nan*ones(md.mesh.numberofvertices, 1);
    md.stressbalance.spcvx=nan*ones(md.mesh.numberofvertices, 1);
    md.stressbalance.spcvy=nan*ones(md.mesh.numberofvertices, 1);
    md.stressbalance.spcvz=nan*ones(md.mesh.numberofvertices, 1);
    md.stressbalance.loadingforce=zeros(md.mesh.numberofvertices, 3);
    md.stressbalance.referential=nan*ones(md.mesh.numberofvertices, 6);
    md.initialization.vx=zeros(md.mesh.numberofvertices,1);
    md.initialization.vy=zeros(md.mesh.numberofvertices,1);
    md.initialization.vz=zeros(md.mesh.numberofvertices,1);
    md.initialization.vel=zeros(md.mesh.numberofvertices,1);
    md.initialization.pressure=md.materials.rho_ice*md.constants.g*md.geometry.thickness;

    % boundary conditions (no-slip at ice divide, free-slip at lateral walls, water pressure at ice front)
    rough_eq=@(x,y) abs(x-y)<10;
    ice_divide=rough_eq(md.mesh.x, min(md.mesh.x));
    ice_front=(md.mesh.vertexonboundary & rough_eq(md.mesh.x, max(md.mesh.x)));
    lateral_wall=md.mesh.vertexonboundary ... 
                 & (rough_eq(md.mesh.y, max(md.mesh.y)) | rough_eq(md.mesh.y, min(md.mesh.y)));
    md.stressbalance.spcvx(ice_divide)=0.0;
    md.stressbalance.spcvy(ice_divide)=0.0;
    md.stressbalance.spcvz(ice_divide)=0.0;
    md.stressbalance.spcvy(lateral_wall)=0.0;
    if NO_SLIP 
        md.stressbalance.spcvx(lateral_wall)=0.0;
    end
    md.mask.ice_levelset(ice_front)=0.0;
    % reinitializaing the ice levelset forces it to be a signed distance function
    if ICE_TONGUE_X<0
        md.mask.ice_levelset=reinitializelevelset(md,md.mask.ice_levelset);
    end
    md.levelset.spclevelset(ice_front)=md.mask.ice_levelset(ice_front);

    if ICE_TONGUE_X>0
        disp(['[init_BC]::Creating an ice tongue starting at x=' num2str(round(ICE_TONGUE_X/1000)) 'km']);
        icetongue_walls = lateral_wall & (md.mesh.x>=ICE_TONGUE_X);
        % md.stressbalance.spcvy(icetongue_walls)=nan;
        % md.stressbalance.spcvx(icetongue_walls)=nan;
        md.mask.ice_levelset(icetongue_walls)=0.;
        md.mask.ice_levelset=reinitializelevelset(md,md.mask.ice_levelset);
        md.levelset.spclevelset(icetongue_walls)=md.mask.ice_levelset(icetongue_walls);
    end
end