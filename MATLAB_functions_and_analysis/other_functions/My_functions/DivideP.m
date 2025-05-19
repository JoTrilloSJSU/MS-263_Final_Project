
   function increasingSegments = DivideP(arr, tolerance, minLength)
    % Initialize variables
    n = length(arr);
    increasingSegments = {};
    segmentStart = 1;
    i = 1;
    
    while i < n
        % Check if the current segment is increasing or within tolerance
        if arr(i) < arr(i + 1) || (arr(i) - arr(i + 1) <= tolerance)
            i = i + 1;
        else
            % Segment has ended
            if (i - segmentStart + 1) >= minLength
                increasingSegments{end + 1} = arr(segmentStart:i);
            end
            % Start a new segment
            segmentStart = i + 1;
            i = segmentStart;
        end
    end
    
    % Check the last segment
    if (i - segmentStart + 1) >= minLength
        increasingSegments{end + 1} = arr(segmentStart:i);
    end
end