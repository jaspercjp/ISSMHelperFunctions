function [time, qg_x_int, qg_x_avg, x_gl_avg, Ls] = calc_gl_flux_t(md, varargin)
    p = inputParser;
    addParameter(p, 'SearchL', 200, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    SEARCH_L=p.Results.SearchL;

    Ls=nan*ones(1,length(md.results.TransientSolution));
    x_gl_avg=nan*ones(1,length(md.results.TransientSolution));
    qg_x_int=nan*ones(1,length(md.results.TransientSolution));
    qg_x_avg=nan*ones(1,length(md.results.TransientSolution));
    time=nan*ones(1,length(md.results.TransientSolution));
    % mean_iceshelf_vx = nan*ones(1,length(md.results.TransientSolution));

    for t=1:length(md.results.TransientSolution)
        time(t)=md.results.TransientSolution(t).time - md.results.TransientSolution(1).time;
        sol=md.results.TransientSolution(t);
        if any(strcmp(fieldnames(sol), 'MeshX'))
            x = sol.MeshX; y=sol.MeshY; 
        else
            x=md.mesh.x; y=md.mesh.y;
        end
        [x_gl, y_gl, q_gl_x, ~] = ISSMGLUtils.calc_gl_flux(x, y, sol, 'SearchL', SEARCH_L);
        
        iceshelf_idx=find(sol.MaskOceanLevelset<0);
        % mean_iceshelf_vx(t)=mean(sol.Vx(iceshelf_idx));
        Ls_tmp=max(x(iceshelf_idx)) - min(x(iceshelf_idx));
        if isempty(Ls_tmp)
            break
        end
        Ls(t)=Ls_tmp;
        x_gl_avg(t)=mean(x_gl);
        qg_x_avg(t)=mean(q_gl_x);
        qg_x_int(t) = trapz(y_gl, q_gl_x);
    end
end