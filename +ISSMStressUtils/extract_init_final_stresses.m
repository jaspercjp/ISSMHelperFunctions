
function [R, R_IT, vel, H, t, last_available_solution_num] = extract_init_final_stresses(md, axis)
    MIN_THICKNESS=-1;
    rho_ice = md.materials.rho_ice; 
    rho_water = md.materials.rho_water;
    g=md.constants.g;
    delta=1-rho_ice/rho_water; 
    % -----------------------------------------
    num_steps = 2;
    t={}; H={};
    Ryy={}; Ryy_IT={};
    sxx={}; syy={};
    vel={};
    idx=0;

    if md.transient.amr_frequency>0
        for i=length(md.results.TransientSolution):-1:1
            if isfield(md.results.TransientSolution(i), 'MeshX')
                last_available_solution_num=i; 
                break;
            end
        end
    else
        last_available_solution_num=length(md.results.TransientSolution);
    end
 
    for i=[1, last_available_solution_num]
        idx=idx+1;
        sol = md.results.TransientSolution(i);
        H_tmp = sol.Thickness;
        vel_tmp = sol.Vel; 
        syy_tmp = sol.DeviatoricStressyy;
        sxx_tmp = sol.DeviatoricStressxx;
        if isempty(syy_tmp) || isempty(sxx_tmp)
            warning("Empty arrays found");
            Ryy_tmp = nan*ones(size(md.mesh.x));
            Rxx_tmp = nan*ones(size(md.mesh.x));
        else
            Ryy_tmp = 2*syy_tmp + sxx_tmp;
            Rxx_tmp = 2*sxx_tmp + syy_tmp;
        end
        if numel(sxx_tmp)~=numel(syy_tmp)
            fprintf("sxx_tmp size=%d | syy_tmp size=%d\n", numel(sxx_tmp), numel(syy_tmp));
        end
        R_IT_tmp = rho_ice*g*(1-rho_ice/rho_water)*H_tmp/2;
        % Calculate average over floating area
        if isfield(sol, 'MaskIceLevelset')
            % invalid_idx=sol.MaskIceLevelset>0 | H_tmp<MIN_THICKNESS | md.mesh.y<=800000 | md.mesh.y>=1000000 | sol.MaskOceanLevelset>=0 | abs(md.mesh.x)>=2000;
            invalid_idx=sol.MaskIceLevelset>0 | H_tmp<MIN_THICKNESS %| sol.MaskOceanLevelset>=0;
            H_tmp(invalid_idx)=nan;
            syy_tmp(invalid_idx)=nan;
            sxx_tmp(invalid_idx)=nan;
            Ryy_tmp(invalid_idx)=nan;
            Rxx_tmp(invalid_idx)=nan;
            R_IT_tmp(invalid_idx)=nan;
            vel_tmp(invalid_idx)=nan;
        end
        if axis==1
            R{idx}=Rxx_tmp;
        elseif axis==2
            R{idx}=Ryy_tmp;
        end
        vel{idx} = vel_tmp;
        R_IT{idx}=R_IT_tmp;
        t{idx} = sol.time - md.results.TransientSolution(1).time;
    end
end