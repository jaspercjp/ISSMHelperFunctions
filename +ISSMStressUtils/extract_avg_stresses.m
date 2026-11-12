function [Ryy, Ryy_confined, Ryy_IT, Ryy_IT_confined, Ryy_norm, H, t] = extract_avg_stresses(md)
% EXTRACT_AVG_STRESSES is used to compute average stresses over floating ice and over time, 
% further separating it into confined and unconfined regions.
    MIN_THICKNESS=-1;
    rho_ice = md.materials.rho_ice; 
    rho_water = md.materials.rho_water;
    g=md.constants.g;
    delta=1-rho_ice/rho_water; 
    % -----------------------------------------
    num_steps = length(md.results.TransientSolution);
    t = nan*ones(1,num_steps);
    H = nan*ones(1,num_steps);
    Ryy = nan*ones(1,num_steps);
    Ryy_IT = nan*ones(1,num_steps);
    Ryy_confined=nan*ones(1,num_steps);
    Ryy_IT_confined=nan*ones(1,num_steps);
    sxx = nan*ones(1,num_steps);
    syy = nan*ones(1,num_steps);
    for i=1:num_steps
        sol = md.results.TransientSolution(i);
        H_tmp = sol.Thickness;
        syy_tmp = sol.DeviatoricStressyy;
        sxx_tmp = sol.DeviatoricStressxx;
        if isempty(syy_tmp) || isempty(sxx_tmp)
            warning("Empty arrays found");
            Ryy_tmp = nan*ones(size(md.mesh.x));
        else
            Ryy_tmp = 2*syy_tmp + sxx_tmp;
        end
        if numel(sxx_tmp)~=numel(syy_tmp)
            fprintf("sxx_tmp size=%d | syy_tmp size=%d\n", numel(sxx_tmp), numel(syy_tmp));
        end
        Ryy_IT_tmp = rho_ice*g*(1-rho_ice/rho_water)*H_tmp/2;
        % Calculate average over floating area
        if isfield(sol, 'MaskIceLevelset')
            % invalid_idx=sol.MaskIceLevelset>0 | H_tmp<MIN_THICKNESS | md.mesh.y<=800000 | md.mesh.y>=1000000 | sol.MaskOceanLevelset>=0 | abs(md.mesh.x)>=2000;
            invalid_idx=sol.MaskIceLevelset>0 | H_tmp<MIN_THICKNESS | sol.MaskOceanLevelset>=0;
            H_tmp(invalid_idx)=nan;
            syy_tmp(invalid_idx)=nan;
            sxx_tmp(invalid_idx)=nan;
            Ryy_tmp(invalid_idx)=nan;
            Ryy_IT_tmp(invalid_idx)=nan;
        end
        if i==num_steps
            disp(nanmin(H_tmp,[],'all'));
        end
        H(i) = nanmean(H_tmp);
        Ryy_confined_tmp = Ryy_tmp;
        Ryy_IT_confined_tmp = Ryy_IT_tmp;
        Ryy_confined_tmp(md.mesh.y>=1000000)=nan;
        Ryy_IT_confined_tmp(md.mesh.y>=1000000)=nan;
        Ryy_norm_tmp = Ryy_tmp./Ryy_IT_tmp;
        Ryy_norm(i) = nanmean(Ryy_norm_tmp);
        Ryy(i) = nanmean(Ryy_tmp);
        Ryy_IT(i) = nanmean(Ryy_IT_tmp);
        Ryy_confined(i) = nanmean(Ryy_confined_tmp);
        sxx(i) = nanmean(sxx_tmp); syy(i) = nanmean(syy_tmp);
        Ryy_IT_confined(i) = nanmean(Ryy_IT_confined_tmp);
        t(i) = sol.time;
    end
    t = t - t(1);
end