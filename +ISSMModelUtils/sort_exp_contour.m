function sort_exp_contour(inputFile, outputFile, direction)
% sort_exp_file Reads an ISSM .exp file, sorts (x,y) pairs, and writes new file
%
% Usage:
%   sort_exp_file('Antarctica.exp', 'Antarctica_sorted.exp', 'ccw');
%
% inputFile  : path to input .exp file
% outputFile : path to output .exp file
% direction  : 'cw' for clockwise, 'ccw' for counterclockwise

    %--- Step 1: Read the exp file ---
    fid = fopen(inputFile, 'r');
    if fid == -1
        error('Could not open input file.');
    end
    rawLines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    rawLines = rawLines{1};

    % Find the line numbers where coordinates begin (skip header)
    coords = [];
    headerLines = {};
    footerLines = {};
    readingCoords = false;

    for i = 1:length(rawLines)
        line = strtrim(rawLines{i});
        if isempty(line)
            continue;
        end

        nums = sscanf(line, '%f %f');
        if numel(nums) == 2
            coords(end+1,:) = nums'; %#ok<AGROW>
            if ~readingCoords
                headerLines = rawLines(1:i-1);
                readingCoords = true;
            end
        else
            if readingCoords
                footerLines = rawLines(i:end);
                break;
            end
        end
    end

    if isempty(coords)
        error('No coordinates found in the exp file.');
    end

    %--- Step 2: Sort the polygon points ---
    % Compute centroid
    cx = mean(coords(:,1));
    cy = mean(coords(:,2));

    % Compute angles of each point w.r.t. centroid
    angles = atan2(coords(:,2) - cy, coords(:,1) - cx);

    % Sort by angle
    [~, order] = sort(angles);
    coordsSorted = coords(order,:);

    % Ensure direction
    areaSigned = polyarea(coordsSorted(:,1), coordsSorted(:,2));
    if strcmpi(direction, 'cw')
        if areaSigned > 0  % polyarea returns positive for ccw
            coordsSorted = flipud(coordsSorted);
        end
    elseif strcmpi(direction, 'ccw')
        if areaSigned < 0
            coordsSorted = flipud(coordsSorted);
        end
    else
        error('direction must be ''cw'' or ''ccw''');
    end

    %--- Step 3: Write to new exp file ---
    fid = fopen(outputFile, 'w');
    if fid == -1
        error('Could not open output file for writing.');
    end

    % Write header
    for i = 1:length(headerLines)
        fprintf(fid, '%s\n', headerLines{i});
    end

    % Write coordinates
    for i = 1:size(coordsSorted,1)
        fprintf(fid, '%.6f %.6f\n', coordsSorted(i,1), coordsSorted(i,2));
    end

    % Write footer if present
    for i = 1:length(footerLines)
        fprintf(fid, '%s\n', footerLines{i});
    end

    fclose(fid);

    fprintf('Sorted exp file written to %s\n', outputFile);
end
