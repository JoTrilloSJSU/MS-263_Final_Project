function cell_result = temparray( location_vector,temperatures)
    % Initialize a 2x7 array filled with NaN (or any placeholder)
  % Get unique values
  x= location_vector;
  y = temperatures;
unique_values = unique(x);

% Initialize a cell array with NaN values and 3 columns
num_rows = length(unique_values);
num_cols = 3;
cell_result = cell(num_rows, num_cols);

% Fill the cell array
for i = 1:num_rows
    indices = find(x == unique_values(i));
    for j = 1:min(length(indices), num_cols)
        cell_result{i, j} = y{unique_values(i)};
    end
end

disp(cell_result);
% function cell_result = temparray(location_vector, temperatures)
%     % Get unique values in location_vector
%     unique_values = unique(location_vector);
% 
%     % Initialize a cell array with NaN values and 3 columns
%     num_rows = length(unique_values);
%     num_cols = 3;
%     cell_result = cell(num_rows, num_cols);
% 
%     % Loop through each unique location
%     for i = 1:num_rows
%         % Find indices for the current unique location
%         indices = find(location_vector == unique_values(i));
% 
%         % Fill the cell array with corresponding temperatures
%         for j = 1:min(length(indices), num_cols)
%             cell_result{i, j} = temperatures{indices(j)};
%         end
% 
%         % Fill remaining cells with NaN if there are fewer temperatures
%         for k = j:num_cols
%             cell_result{i, k} = NaN; % or some placeholder if needed
%         end
%     end
% 
%     % Display the result
%     disp(cell_result);
% end
