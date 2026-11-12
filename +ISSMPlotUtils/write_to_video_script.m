function path = write_to_video_script()
    % Path of the frames-to-video helper the animate_* functions shell out to.
    % It is not part of this library; point $ISSM_WRITE_TO_VIDEO at your copy.
    path = getenv('ISSM_WRITE_TO_VIDEO');
    if isempty(path)
        error('ISSMHelperFunctions:noVideoScript', ...
              'Set ISSM_WRITE_TO_VIDEO to the write_to_video.py script used to encode the frames.');
    end
end
