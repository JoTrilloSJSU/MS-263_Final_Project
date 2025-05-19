function plot = Ptemp(P, JAC_T, setupfilestr, fs_fast, a)
    % PTEMP creates a plot of pressure and temperature over time (seconds) 
    % where temperature is represented by a color bar. 
    % a = start time (s) (start with default a=0)
    % Must have saved Data_xxxx as .mat file saved using odas_p2mat

    % Find rate of samples taken in samples/sec
    rate = channel_sampling_rate('P', setupfilestr, fs_fast);
    rateT = channel_sampling_rate('JAC_T', setupfilestr, fs_fast);

    % Check if rate of P and JAC_T (temperature) are the same
    if rate ~= rateT
        fprintf('Pressure and Temperature have different sampling rates\n');
        return; % Exit if sampling rates are not the same
    end

    % Total number of samples
    totalsamples = length(P);

    % Time vector
    time = (0:totalsamples - 1) / rate;

    % Adjust start time if provided
    if a ~= 0
        start_index = round(a * rate) + 1; % +1 to handle MATLAB 1-based indexing
        if start_index > totalsamples
            fprintf('Start time exceeds the duration of the data\n');
            return;
        end
        P = P(start_index:end);
        JAC_T = JAC_T(start_index:end);
        time = time(start_index:end);
    end

    % Check if length of P and JAC_T are the same
    if length(JAC_T) ~= length(P)
        fprintf('Pressure and Temperature have different lengths\n');
        return;
    end

    % Create a logical index to filter out P < 1
    valid_indices = P >= 1;

    % Apply the logical index to filter data
    P_filtered = P(valid_indices);
    JAC_T_filtered = JAC_T(valid_indices);
    time_filtered = time(valid_indices);

    % Scatter plot
    scatter(time_filtered, P_filtered, [], JAC_T_filtered, 'fill');
    xlabel('Time (seconds)');
    ylabel('Pressure (dbar)');
    set(gca, 'YDir', 'reverse');

    % Add colorbar
    cb = colorbar;
    colormap(turbo); % Use the turbo colormap
    clim([min(JAC_T_filtered) 14]);

    % Set colorbar label
    ylabel(cb, 'Temperature ({\circ}C)');
end
