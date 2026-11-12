function write_exp(filename, coords, contourName)
% write_exp Write (x,y) coordinates to an ISSM .exp file
%
% Usage:
%   write_exp('myoutline.exp', coords, 'MyContour');
%
% filename    : name of output .exp file
% coords      : Nx2 array of [x y] coordinates
% contourName : string name for the contour (optional)

    if nargin < 3
        contourName = 'DefaultContour';
    end

    % Ensure coords is closed (last point = first point)
    if ~isequal(coords(1,:), coords(end,:))
        coords(end+1,:) = coords(1,:);
    end

    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open %s for writing.', filename);
    end
    
    % Write headers
    fprintf(fid, '## Name:%s\n', contourName);
    fprintf(fid, '## Icon:0\n');
    fprintf(fid, '# Points Count Value\n', size(coords,1));
    fprintf(fid, '%d 1.0\n', size(coords,1));
    fprintf(fid, '# X pos Y pos\n', size(coords,1));

    % Write coordinates
    for i = 1:size(coords,1)
        fprintf(fid, '%.6f %.6f\n', coords(i,1), coords(i,2));
    end

    % Footer
    fprintf(fid, '## End\n');

    fclose(fid);
    fprintf('Wrote %d points to %s\n', size(coords,1), filename);
end
