function temperatures = tempcomparison(P, JAC_T, setupfilestr, fs_fast, a)
    % TEMP_COMPARISON returns the temperature at specified depths of interest
    % a = start time (s) (default a = 0)
    % Must have saved Data_xxxx as .mat file saved using odas_p2mat

    % Define specific depths of interest
    depths_of_interest = 1:40; % Depths from 1 to 40 dbar

    % Find rate of samples taken in samples/sec
    rate = channel_sampling_rate('P', setupfilestr, fs_fast);
    rateT = channel_sampling_rate('JAC_T', setupfilestr, fs_fast);

    % Check if rate of P and JAC_T (temperature) are the same
    if rate ~= rateT
        fprintf('Pressure and Temperature have different sampling rates\n');
        temperatures = [];
        return; % Exit if sampling rates are not the same
    end

    % Total number of samples
    totalsamples = length(P);

    % Adjust start time if provided
    if a ~= 0
        start_index = round(a * rate) + 1; % +1 to handle MATLAB 1-based indexing
        if start_index > totalsamples
            fprintf('Start time exceeds the duration of the data\n');
            temperatures = [];
            return;
        end
        P = P(start_index:end);
        JAC_T = JAC_T(start_index:end);
    end

    % Check if length of P and JAC_T are the same
    if length(JAC_T) ~= length(P)
        fprintf('Pressure and Temperature have different lengths\n');
        temperatures = [];
        return;
    end

    % Initialize the temperatures array
    temperatures = nan(length(depths_of_interest), 1);

    % Print the range of depths available
    fprintf('Depth range in data: %.2f to %.2f dbar\n', min(P), max(P));

    % Find temperatures at closest depths
    for i = 1:length(depths_of_interest)
        depth = depths_of_interest(i);
        
        % Check if the depth is within the range of P
        if depth < min(P) || depth > max(P)
            fprintf('Depth %d dbar is outside the range of data\n', depth);
            continue;
        end
        
        % Find the index of the closest depth
        [~, idx] = min(abs(P - depth)); % Find index of closest value
        
        if ~isempty(idx)
            % Store temperature at the closest depth
            temperatures(i) = JAC_T(idx);
        else
            fprintf('No data found for depth %d dbar\n', depth);
        end
    end

    % Display the results
    for i = 1:length(depths_of_interest)
        if ~isnan(temperatures(i))
            fprintf('Depth: %d dbar, Temperature: %.2f°C\n', depths_of_interest(i), temperatures(i));
        else
            fprintf('Depth: %d dbar, Temperature: Data not available\n', depths_of_interest(i));
        end
    end
end