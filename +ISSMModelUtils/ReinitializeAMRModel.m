%REINITIALIZEAMRMODEL Reinitializes an ISSM model after AMR transient simulation.
%
%   md = REINITIALIZEAMRMODEL(md, config)
%   md = REINITIALIZEAMRMODEL(md, config, 'ParamFunc', @(md) md)
%
%   This function updates the ISSM model structure after a transient simulation,
%   especially when adaptive mesh refinement (AMR) has been used. It copies the
%   final mesh and solution fields, optionally applies a parametrization function,
%   and resets initial conditions and geometry for further simulations.
%
%   Inputs:
%       md         - ISSM model structure containing transient results.
%       config     - Configuration structure (currently unused).
%       ParamFunc  - (Optional) Function handle to further parametrize the model
%                    after mesh update (default: @(x) nan).
%
%   Outputs:
%       md         - ISSM model structure reinitialized with updated mesh,
%                    geometry, and initial fields.
%
%   The function also moves previous transient results to a new field in md.results
%   and resets the start time for the next simulation.
%
%   Example:
%       md = ISSMModelUtils.ReinitializeAMRModel(md, [], 'ParamFunc', @my_param_func)
function md = ReinitializeAMRModel(md, varargin)
    p = inputParser;
    addParameter(p, 'ParamFunc', @(x) nan, @(x) isa(x, 'function_handle'));
    addParameter(p, 'IceCutCD', nan, @(x) isnan(x)||isnumeric(x));
    addParameter(p, 'TimeIdx', nan, @(x) isnan(x)||isnumeric(x) )
    parse(p,varargin{:});
    
    PARAM_FUNC=p.Results.ParamFunc;
    ICE_CUT_CD=p.Results.IceCutCD;
    TIME_IDX=p.Results.TimeIdx;

    assert(isfield(md.results,'TransientSolution'), 'ERROR: Input to ReinitializeModel must be a transient simulation.');

    if ~isnan(TIME_IDX)
        sol = md.results.TransientSolution(TIME_IDX);
    else
        sol = md.results.TransientSolution(end);
    end

    if isempty(sol)
        sol=md.results.TransientSolution(end-1);
    elseif isfield(sol,'MeshX')
        if isempty(sol.MeshX)
            sol=md.results.TransientSolution(end-1);
        end
    end
    
    % Check if amr was used to run the model 
    if isfield(md.results.TransientSolution(1),'MeshX')
        isamr= true;
    else
        isamr=false;
    end
    
    % No need to re-initialize the mesh if AMR was not used.
    if isamr
        disp('- [ReinitializeAMRModel] :: copying amr mesh...')
        md=ISSMMeshUtils.copy_amr_mesh(md, sol);
        md.mesh
        md.geometry
        md_tmp=PARAM_FUNC(md); % Call PARAM_FUNC to evaluate the geometry on the new grid
        if isa(md_tmp, 'model')
            md=md_tmp;
        end
    else
        md_tmp=PARAM_FUNC(md);
        if ~isa(md_tmp, 'model') && ~isnan(md_tmp)
            warning('- [ReinitializeAMRModel] :: Improper parametrization function was given');
        end
    end
    
    disp('- [ReinitializeAMRModel] :: copying geometry...')
    if isfield(sol, 'Bed')
        md.geometry.surface = sol.Surface; 
        md.geometry.base = sol.Base;
        md.geometry.thickness = sol.Thickness;
        md.geometry.bed = sol.Bed;
        md.mask.ocean_levelset = sol.MaskOceanLevelset;
    else
        md.geometry.thickness = sol.Thickness;
        md.geometry.base = md.geometry.bed;
        floating_base=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness;
        pos = find(floating_base>md.geometry.bed);
        md.geometry.base(pos) = floating_base(pos);
        md.geometry.surface = md.geometry.base + md.geometry.thickness;
        di=md.materials.rho_ice/md.materials.rho_water;
        md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;
    end
    
    md.initialization.vx = sol.Vx;
    md.initialization.vy = sol.Vy; 
    md.initialization.vel = sol.Vel;
    md.initialization.pressure=sol.Pressure;
    % md.results.reinitialization.DeviatoricStressxx = sol.DeviatoricStressxx;
    % md.results.reinitialization.DeviatoricStressyy = sol.DeviatoricStressyy; 
    md.results.reinitialization.Thickness = sol.Thickness;
    md.initialization.pressure=sol.Pressure;
    % md.mask.ice_levelset=sol.MaskIceLevelset;

    % Modify ice levelset if a calving criterion is provided
    if ~isnan(ICE_CUT_CD)
        rho_ice=md.materials.rho_ice; rho_water=md.materials.rho_water;
        g=md.constants.g;
        Rxx=2*sol.DeviatoricStressxx + sol.DeviatoricStressyy;
        Rxx_IT=rho_ice*(1-rho_ice/rho_water)*g*sol.Thickness/2;
        Rxx_norm = Rxx ./ Rxx_IT; K=1-Rxx_norm;
        cd=1-sqrt(max(K, 0));
        to_cut = ISSMModelUtils.find_contiguous_calving_nodes(md, cd>=ICE_CUT_CD);
        md.mask.ice_levelset(to_cut & md.mask.ocean_levelset<0)=1;
        % md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);
    end

    md=setflowequation(md,'SSA','all');
    md.timestepping.start_time = md.timestepping.final_time;

    % Allocate a new field in results for the outputs of the new simulation
    newname = ['TransientSolution' num2str(numel(fields(md.results))+1)];
    if isfield(md.results,newname)
        error(['- [ReinitializeAMRModel] :: Cannot save ' newname ' in md.results']);
    else
        disp(['- [ReinitializeAMRModel] :: Moving results to ' newname]);
        md.results.(newname) = md.results.TransientSolution;
        md.results.TransientSolution  = struct();
    end
end
