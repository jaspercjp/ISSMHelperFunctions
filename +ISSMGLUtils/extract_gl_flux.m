function [q, v, H, t] = extract_gl_flux(md, axis)
    MIN_THICKNESS=-1;
    % -----------------------------------------
    num_steps = length(md.results.TransientSolution);
    t = nan*ones(1,num_steps);
    q = nan*ones(1,num_steps);
    v = nan*ones(1,num_steps);
    H = nan*ones(1,num_steps);
    
    for i=1:num_steps
        sol = md.results.TransientSolution(i);
        H_tmp = sol.Thickness;
        Vy_tmp = sol.Vy;
        Vx_tmp = sol.Vx; 
        
        % % Calculate average over floating area
        % if isfield(sol, 'MaskIceLevelset')
        %     % invalid_idx=sol.MaskIceLevelset>0 | H_tmp<MIN_THICKNESS | md.mesh.y<=800000 | md.mesh.y>=1000000 | sol.MaskOceanLevelset>=0 | abs(md.mesh.x)>=2000;
        %     invalid_idx=sol.MaskIceLevelset>0 | H_tmp<MIN_THICKNESS;
        %     H_tmp(invalid_idx)=nan;
        %     Vy_tmp(invalid_idx)=nan;
        % end
        
        if axis==1
            v_long = Vx_tmp; 
        elseif axis==2
            v_long = Vy_tmp;
        end
        q_all = v_long .* H_tmp; 
        gl_idx=abs(sol.MaskOceanLevelset) < 300;
        if all(sol.MaskOceanLevelset(gl_idx)<0)
            disp(t(i-1))
            warning('No Grounded Ice Found in Model');
            break;
        end
        q_gl = q_all(gl_idx);
        num_gl_nodes = numel(q_gl);
        q(i) = nanmean(q_gl);
        v(i) = nanmean(v_long(gl_idx));
        H(i) = nanmean(H_tmp(gl_idx));
        t(i) = sol.time;
    end
    t = t - t(1);
end