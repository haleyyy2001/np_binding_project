% ... [previous code to generate data_table] ...
% Define levels for each parameter

%tcr_params = [1000, 15, 50, 20, 300]; % tcr_params_cluster
random_numbers = rand(4, 1);
random_numsmall = rand(3, 1);
% Scale the numbers to the range [min_value, max_value] 
%tcr_params = [1000, 15, 50, 20, 300]; 
 tcr_params = [1000, 300, 50, 1, 300]; 

%vh_levels = 10+ 90 * random_numbers; % linspace(10, 100, 15);


vh_levels=[1,2,3,4,5,10,15,20,25,30,50,70,100];
np_radius_levels = 10+ 70 * random_numbers;%linspace(10, 80, 15);
a = -1+ 3 * random_numsmall  ;  %linspace(-1, 2, 15);
koff_levels = 10.^a;

% Factorial design
[DOE_vh, DOE_np_radius, DOE_koff] = ndgrid(vh_levels, np_radius_levels, koff_levels);
full_experiments = numel(DOE_vh);

% Selecting 30% of the points randomly
selected_indices = randperm(full_experiments, round(full_experiments * 1));

% Reshape the DOE arrays into vectors
DOE_vh_vector = DOE_vh(:);
DOE_np_radius_vector = DOE_np_radius(:);
DOE_koff_vector = DOE_koff(:);

% Initialize arrays for the selected experiments using reshaped vectors
selected_DOE_vh = DOE_vh_vector(selected_indices);
selected_DOE_np_radius = DOE_np_radius_vector(selected_indices);
selected_DOE_koff = DOE_koff_vector(selected_indices);
DOE_bind_freq = zeros(length(selected_indices), 1);

% Run experiments for the selected indices
for idx = 1:length(selected_indices)
    i = selected_indices(idx);
    np_params = [selected_DOE_vh(idx), selected_DOE_np_radius(idx), 1];
    kinetic_params = [0.01, 0.1, selected_DOE_koff(idx)];
    [bind_freq, nt_bound, nt_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0);
    DOE_bind_freq(idx) = sum(bind_freq);
end

% Preparing data for parallelplot
data_table = table(selected_DOE_vh, selected_DOE_np_radius, selected_DOE_koff, DOE_bind_freq, ...
                   'VariableNames', {'vh', 'np radius', 'koff', 'bindprobability'});
% ... [previous code to generate data_table] ...

% Normalization function
normalize = @(x) (x - min(x)) / (max(x) - min(x));

% Apply normalization to each column of the data table
normalized_data_table = varfun(normalize, data_table);

% Normalize bind_probability for color mapping
normalized_bind_prob = normalize(data_table.bindprobability);

% Define the color mapping function
getColor = @(x) [1, 0.6, 0.8] * (x < 0.5) + ...
                [0.6/2, 0.85/2, 1/2] * (x >= 0.5 & x <= 1) ;
               % [0.6, 0.4, 0.8] * (x >= 0.7);

% Convert the table to an array for plotting
data_array = table2array(normalized_data_table);

% Create a new figure for the parallel coordinates plot
figure;
hold on; % Enable adding multiple lines to the plot

% Plot each line individually with its color
for i = 1:size(data_array, 1)
    plot(1:size(data_array, 2), data_array(i, :), 'Color', getColor(normalized_bind_prob(i)));
end

% Set the labels for the axes
set(gca, 'xtick', 1:size(data_array, 2), 'xticklabel',[  {'vh'}    {' np radius'}    {'koff'}    {'bind probability'} ]);



% Find the y-axis limits for normalization
y_limits = get(gca, 'YLim');

% Draw a vertical line at the position of np_radius
koff_position = 1;
line([koff_position, koff_position], y_limits, 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');

% Add scale on the line and draw horizontal grid lines
scale_values =  (data_table.vh);
normalized_scale_values = (scale_values - min(scale_values)) / (max(scale_values) - min(scale_values)); % Normalize scale values to [0, 1]
normalized_scale_values = normalized_scale_values * (y_limits(2) - y_limits(1)) + y_limits(1); % Map to actual y-limits
 
for i = 1:length(scale_values)
    % Add text labels
    text(koff_position + 0.03, normalized_scale_values(i), num2str((scale_values(i))), 'HorizontalAlignment', 'left');
    
    % Draw horizontal grid lines across the np_radius bar
    line([koff_position- 0.03, koff_position + 0.03], [normalized_scale_values(i), normalized_scale_values(i)], 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%


% Find the y-axis limits for normalization
y_limits = get(gca, 'YLim');

% Draw a vertical line at the position of np_radius
np_radius_position = 2;
line([np_radius_position, np_radius_position], y_limits, 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');

% Add scale on the line and draw horizontal grid lines
scale_values = 10:10:80; % Scale from 10 to 80 with a gap of 10
normalized_scale_values = (scale_values - min(scale_values)) / (max(scale_values) - min(scale_values)); % Normalize scale values to [0, 1]
normalized_scale_values = normalized_scale_values * (y_limits(2) - y_limits(1)) + y_limits(1); % Map to actual y-limits

for i = 1:length(scale_values)
    % Add text labels
    text(np_radius_position + 0.03, normalized_scale_values(i), num2str(scale_values(i)), 'HorizontalAlignment', 'left');
    
    % Draw horizontal grid lines across the np_radius bar
    line([np_radius_position - 0.03, np_radius_position + 0.03], [normalized_scale_values(i), normalized_scale_values(i)], 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');
end

% ... [rest of your code] ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%%%
% ... [previous plotting code] ...




% Find the y-axis limits for normalization
y_limits = get(gca, 'YLim');

% Draw a vertical line at the position of koff
koff_position = 3;
line([koff_position, koff_position], y_limits, 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');

% Add scale on the line and draw horizontal grid lines
scale_values = -1:0.25:2; % Range of the log-transformed variable 'a'

% Calculate the koff_levels from scale_values
koff_levels = 10.^scale_values;

% Normalize koff_levels to [0, 1]
normalized_koff_levels = (koff_levels - min(koff_levels)) / (max(koff_levels) - min(koff_levels));

% Map normalized_koff_levels to actual y-limits
mapped_koff_levels = normalized_koff_levels * (y_limits(2) - y_limits(1)) + y_limits(1);

for i = 1:length(koff_levels)
    % Add text labels for koff_levels
    text(koff_position + 0.03, mapped_koff_levels(i), num2str(koff_levels(i)), 'HorizontalAlignment', 'left');
    
    % Draw horizontal grid lines across the koff bar
    line([koff_position - 0.03, koff_position + 0.03], [mapped_koff_levels(i), mapped_koff_levels(i)], 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');
end













 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%%%
% ... [previous plotting code] ...

% Find the y-axis limits for normalization
y_limits = get(gca, 'YLim');

% Draw a vertical line at the position of np_radius
np_radius_position = 4;
line([np_radius_position, np_radius_position], y_limits, 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');

% Add scale on the line and draw horizontal grid lines
scale_values = min(data_table.bindprobability):(max(data_table.bindprobability)- min(data_table.bindprobability))/10:max(data_table.bindprobability); % Scale from 10 to 80 with a gap of 10


normalized_scale_values = (scale_values - min(scale_values)) / (max(scale_values) - min(scale_values)); % Normalize scale values to [0, 1]
normalized_scale_values = normalized_scale_values * (y_limits(2) - y_limits(1)) + y_limits(1); % Map to actual y-limits

for i = 1:length(scale_values)
    % Add text labels
    text(np_radius_position + 0.03, normalized_scale_values(i), num2str(scale_values(i)), 'HorizontalAlignment', 'left');
    
    % Draw horizontal grid lines across the np_radius bar
    line([np_radius_position - 0.03, np_radius_position + 0.03], [normalized_scale_values(i), normalized_scale_values(i)], 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');
end

 title("Uniform Surface")
hold off; % Disable adding more lines to the plot











