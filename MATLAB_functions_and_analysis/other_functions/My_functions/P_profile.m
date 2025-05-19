function [profile_pressure] = P_profile(fname, P_start, P_end, varargin)
% P_PROFILE Extracts and returns pressure data for a specific profile
%   fname      : Filename for the data
%   P_start    : Starting pressure for the profile
%   P_end      : Ending pressure for the profile
%   varargin   : Additional parameters for customization

%-----------------------------------------------------------------
% ----- Default parameters ---------------------------------------
%-----------------------------------------------------------------
default_YD_0 = 0;                % Deployment day [year day]
default_make_figures = true;     % Flag to render figures

% -- Values used to determine profiles
default_profile_num = 1;         % Profile number
default_profile_min_P = 1;       % Minimum pressure [dBar]
default_profile_min_W = 0.2;     % Minimum speed [m/s]
default_profile_min_duration = 20; % Minimum duration [s]

% Return defaults if no inputs are specified
if nargin == 0
    for d = whos('default_*')'
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    return
end

%-----------------------------------------------------------------
% ----- Parse Inputs ---------------------------------------------
%-----------------------------------------------------------------
p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'fname', @ischar);
addRequired(p, 'P_start', @isnumeric);
addRequired(p, 'P_end', @isnumeric);

addParameter(p, 'YD_0', default_YD_0, @isnumeric);
addParameter(p, 'make_figures', default_make_figures, @islogical);
addParameter(p, 'profile_num', default_profile_num, @isnumeric);
addParameter(p, 'profile_min_P', default_profile_min_P, @isnumeric);
addParameter(p, 'profile_min_W', default_profile_min_W, @isnumeric);
addParameter(p, 'profile_min_duration', default_profile_min_duration, @isnumeric);

parse(p, fname, P_start, P_end, varargin{:});

% Extract parameters
YD_0 = p.Results.YD_0;
make_figures = p.Results.make_figures;
profile_num = p.Results.profile_num;
profile_min_P = p.Results.profile_min_P;
profile_min_W = p.Results.profile_min_W;
profile_min_duration = p.Results.profile_min_duration;



%-----------------------------------------------------------------
% ----- Load Data ------------------------------------------------
%-----------------------------------------------------------------
d = odas_p2mat(fname, p.Results);
[P, N, E] = fileparts(d.fullPath);
File_Name = [N E];
if ~exist([P filesep N '.mat'], 'file')
    disp(['Saving into MAT-file: ' P filesep N '.mat']);
    save([P filesep N '.mat'], '-struct', 'd', '-v6');
else
    disp(['Loading from MAT-file: ' P filesep N '.mat']);
end

% Extract fields from data structure
for field = fieldnames(d)'
    eval([char(field) ' = d.' char(field) ';']);
end

% Define profile_dir based on the data or default to empty
profile_dir = '';
if exist('cfgobj', 'var')
    profile_dir = char(setupstr(cfgobj, 'instrument_info', 'profile_dir'));
end
if isempty(profile_dir)
    profile_dir = 'up'; % Default value if profile_dir is not available
end


%-----------------------------------------------------------------
% ----- Determine Profile ----------------------------------------
%-----------------------------------------------------------------
profile = [1; length(t_slow)]; % start using entire file

if exist('P_slow', 'var') && exist('W_slow', 'var')
    if strcmpi(profile_dir, 'up') || strcmpi(profile_dir, 'down')
        profile = get_profile(P_slow, W_slow, profile_min_P, ...
            profile_min_W, profile_dir, ...
            profile_min_duration, fs_slow);
    elseif strcmpi(profile_dir, 'glide')
        profile_down = get_profile(P_slow, W_slow, profile_min_P, ...
            profile_min_W, 'down', ...
            profile_min_duration, fs_slow);
        profile_up = get_profile(P_slow, W_slow, profile_min_P, ...
            profile_min_W, 'up', ...
            profile_min_duration, fs_slow);
        % Sort columns in ascending order
        profile = sort([profile_down profile_up], 2);
    end
    
    % Debugging output
    disp(['Size of profile array: ', num2str(size(profile))]);
    disp(['profile_num: ', num2str(profile_num)]);

    % Validate profile_num
    if profile_num > size(profile, 2)
        error('profile_num exceeds the number of columns in the profile array.');
    end


    % Extract pressure data for the specific profile
    start_index_slow = profile(1, profile_num);
    end_index_slow = profile(2, profile_num);
    
    start_index_fast = 1 + round((fs_fast / fs_slow) * (start_index_slow - 1));
    end_index_fast = round((fs_fast / fs_slow) * end_index_slow);
    m = start_index_slow:end_index_slow;
    
    profile_pressure = P_slow(m);
end

end
