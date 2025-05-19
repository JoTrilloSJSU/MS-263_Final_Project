function temperatures = tempcompcells(P_cells, JAC_T_cells, setupfilestr, fs_fast,binsize)
    % TEMP_COMPARISON returns the average temperature at specified depths of interest
    % for multiple sets of pressure and temperature data.
    
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
end
