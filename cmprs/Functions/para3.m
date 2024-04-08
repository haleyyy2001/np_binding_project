% Define levels for each parameter
tcr_params = [1000, 300, 50, 1, 300]; 
vh_levels = linspace(10, 100, 15);
np_radius_levels = linspace(10, 80, 15);
a = linspace(-1, 2, 15);
koff_levels = 10.^a;
 

% Factorial design
[DOE_vh, DOE_np_radius, DOE_koff] = ndgrid(vh_levels, np_radius_levels, koff_levels);
full_experiments = numel(DOE_vh);

% Selecting 30% of the points randomly
selected_indices = randperm(full_experiments, round(full_experiments * 0.3));

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
                   'VariableNames', {'vh', 'np_radius', 'koff', 'bind_probability'});

% Create a parallel coordinates plot with data normalization
figure('Units','normalized','Position',[0.3 0.3 0.45 0.4]);
p = parallelplot(data_table);

% Apply custom colors to the lines after rendering
drawnow; % Ensure the plot is rendered before modifying

% Define the color mapping function
getColor = @(x) [1, 0.6, 0.8] * (x < 0.3) + ...
                [0.68, 0.85, 0.9] * (x >= 0.3 & x < 0.7) + ...
                [0.6, 0.4, 0.8] * (x >= 0.7);

% Calculate colors for each line
lineColors = arrayfun(getColor, data_table.bind_probability, 'UniformOutput', false);

% Apply colors to each line
for i = 1:length(p.NodeChildren)
    child = p.NodeChildren(i);
    if isa(child, 'matlab.graphics.chart.primitive.Line')
        % Directly use the loop index to set the color
        if i <= length(lineColors)
            child.Color = lineColors{i};
        end
    end
end

% Customize the plot if needed
% [Additional customization code here, if necessary]
