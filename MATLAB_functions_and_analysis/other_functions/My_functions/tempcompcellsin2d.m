function temperatures = tempcompcellsin2d(P_cells, JAC_T_cells, setupfilestr, fs_fast,binsize,location_vector)
    % TEMP_COMPARISON returns the average temperature at specified depths of interest
    % for multiple sets of pressure and temperature data.
    
    % Define specific depths of interest
    depths_of_interest = 1:binsize:40; % Depths from 1 to 40 dbar
    
    % Find rate of samples taken in samples/sec
    rate = channel_sampling_rate('P', setupfilestr, fs_fast);
    rateT = channel_sampling_rate('JAC_T', setupfilestr, fs_fast);
    
    % Check if rate of P and JAC_T (temperature) are the same
    if rate ~= rateT
        fprintf('Pressure and Temperature have different sampling rates\n');
        temperatures = [];
        return; % Exit if sampling rates are not the same
    end
    
    % Initialize cell array to hold temperature results for each dataset
    num_datasets = length(P_cells);
    temperatures = cell(num_datasets, 1);
    
    % Process each dataset
    for dataset_idx = 1:num_datasets
        P = P_cells{dataset_idx};
        JAC_T = JAC_T_cells{dataset_idx};
      
        
        % Check if length of P and JAC_T are the same
        if length(JAC_T) ~= length(P)
            fprintf('Pressure and Temperature have different lengths for dataset %d\n', dataset_idx);
            temperatures{dataset_idx} = [];
            continue;
        end
        
        % Initialize array to hold temperature sums and counts
        temp_sums = zeros(length(depths_of_interest), 1);
        temp_counts = zeros(length(depths_of_interest), 1);
        
        % Print the range of depths available
        fprintf('Dataset %d: Depth range in data: %.2f to %.2f dbar\n', dataset_idx, min(P), max(P));
        
        % Find average temperatures at specified depths
        for sample_idx = 1:length(P)
            current_depth = P(sample_idx);
            current_temp = JAC_T(sample_idx);
            
            % Check if current depth is in the specified depths_of_interest
            for i = 1:length(depths_of_interest)
                if current_depth >= depths_of_interest(i) && current_depth < depths_of_interest(i) + 1
                    temp_sums(i) = temp_sums(i) + current_temp; % Accumulate temperature
                    temp_counts(i) = temp_counts(i) + 1; % Increment count
                end
            end
        end
        
        % Compute average temperatures
        dataset_temperatures = nan(length(depths_of_interest), 1);
        for i = 1:length(depths_of_interest)
            if temp_counts(i) > 0
                dataset_temperatures(i) = temp_sums(i) / temp_counts(i); % Calculate average
            end
        end
        
        % Store results in the cell array
        temperatures{dataset_idx} = dataset_temperatures;
        
        % Display the results for this dataset
        for i = 1:length(depths_of_interest)
            if ~isnan(dataset_temperatures(i))
                fprintf('Dataset %d: Depth: %d dbar, Average Temperature: %.2f°C\n', dataset_idx, depths_of_interest(i), dataset_temperatures(i));
            else
                fprintf('Dataset %d: Depth: %d dbar, Average Temperature: Data not available\n', dataset_idx, depths_of_interest(i));
            end
        end
    end

    % Create a 2d array that sorts the location of each profile into
    % columns

  % Find dimensions of array
% unique_locations = unique(location_vector);
% length_temp_array = length(unique_locations); % Get the number of unique locations
% mode_location_array = mode(location_vector);
% height_temp_array = 3; % This should return a scalar

% Parameters
height_temp_array = 3;  % Number of rows
length_temp_array = 7;   % Number of columns
temperatures = cell(14, 1); % Example initialization of temperatures
for i = 1:14
    temperatures{i} = (1:40) * i; % Populate each cell with a 40×1 double
end

% Initialize result array with NaNs
result_array = cell(height_temp_array, length_temp_array);
for row = 1:height_temp_array
    for col = 1:length_temp_array
        result_array{row, col} = nan; % Fill with NaN
    end
end

% Fill the result array based on the described pattern
instance_counter = 1; % Counter for the instance of each number
for value_idx = 1:length(temperatures)
    temperature_data = temperatures{value_idx}; % Get the 40x1 double array
    
    % Loop through the elements of the current temperature data
    for temp_idx = 1:length(temperature_data)
        % Determine the row and column based on the instance
        row = ceil(instance_counter / length_temp_array);
        col = mod(instance_counter - 1, length_temp_array) + 1;

        % Assign the value to the appropriate position
        result_array{row, col} = temperature_data(temp_idx); 
        instance_counter = instance_counter + 1; % Increment the instance counter

        % Stop if we have filled the entire result array
        if instance_counter > (height_temp_array * length_temp_array)
            break;
        end
    end
    if instance_counter > (height_temp_array * length_temp_array)
        break;
    end
end

% Convert to numeric array for better display
result_array_numeric = cell2mat(result_array);
disp(result_array_numeric);

end
