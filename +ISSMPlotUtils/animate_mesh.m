function animate_mesh(md, DATA_DIR ,varargin)
    p = inputParser;

    % model setup
    addParameter(p, 'NumYrs', 200, @(x) isnumeric(x) && x>0);
    addParameter(p, 'Skip', 5, @(x) isnumeric(x) && x>0);
    addParameter(p, 'FrameRate', 20, @(x) isnumeric(x) && x>0);
    parse(p, varargin{:});
    NUM_YRS=p.Results.NumYrs
    SKIP=p.Results.Skip;
    FRAME_RATE=p.Results.FrameRate;

    FRAME_DIR=DATA_DIR+"mesh_frames/";
    if exist(FRAME_DIR,'dir')
        system(sprintf("rm %s*.png", FRAME_DIR));
    end
    if ~exist(FRAME_DIR,'dir')
        system(sprintf("mkdir %s", FRAME_DIR));
    end

    cm_pert = cm_interp([0 0 1; 0.9 0.9 0.9; 1 0 0], 21);
    fig=figure; 
    % tiledlayout(1,1,'TileSpacing','compact','Padding','loose');

    
    num_t = length(md.results.TransientSolution);
    for t=1:SKIP:num_t
        sol = md.results.TransientSolution(t);
        t_pert = sol.time - md.results.TransientSolution(1).time;
        if (t_pert>NUM_YRS)
            break
        end

        % vx1(sol.MaskIceLevelset>0)=nan;
        if isfield(sol, 'MeshX')
            if isempty(sol.MeshX)
                continue;
            end
            set(fig, 'Position', [100, 100, 1700, 550]); % [left, bottom, width, height]
            plotmodel(md,'data','mesh','amr',t)
        else
            error('Model mesh is not changing in time.')
            return;
        end
        

        sgt=sgtitle(sprintf("\\bf Time=%.1f yr", t_pert), ...
            'interpreter','latex', 'FontSize',20, 'Color', 'red');

        % for i=1:4
        %     nexttile(i); xlim([0 320]);
        %     % xline(640, 'k--', 'LineWidth', 3)
        % end
         % Write to video 
        set(gcf, 'Color', 'w');            % Set figure background to white
        
        drawnow; 
        frame = getframe(gcf); % Capture current figure as a frame
        if t<100
            if t<10 
                frame_num="00"+string(t);
            else
                frame_num = "0" + string(t);
            end
        else
            frame_num = string(t);
        end
        imwrite(frame.cdata, FRAME_DIR+sprintf("frame%s.png",frame_num));  % Save as PNG image
    end
    command = sprintf("python3 %s %s %s %d", ISSMPlotUtils.write_to_video_script(), ...
                            FRAME_DIR, DATA_DIR + "mesh_animation.mp4", FRAME_RATE);
    system(command);
end