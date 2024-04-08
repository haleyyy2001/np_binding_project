% Define levels for each parameter
vh_levels = [10, 20, 30];
np_radius_levels = linspace(50, 500, 2);
a = linspace(-1, 2, 3);
koff_levels = 10.^a;

% Factorial design using ndgrid
[DOE_vh, DOE_np_radius, DOE_koff] = ndgrid(vh_levels, np_radius_levels, koff_levels);
num_experiments = numel(DOE_vh);

% Initialize the binding frequency array
DOE_bind_freq = zeros(num_experiments, 1);

% Define tcr_params
tcr_params = [1000, 15, 50, 20, 300];

% Calculate binding frequency for each experiment
for i = 1:num_experiments
    np_params = [DOE_vh(i), DOE_np_radius(i), 1];
    kinetic_params = [0.01, 0.1, DOE_koff(i)];
    bind_freq = NP_surface_binding(np_params, kinetic_params, tcr_params, 0);
    DOE_bind_freq(i) = mean(bind_freq(:));
end

% Normalize only binding frequency data
norm_bind_freq = (DOE_bind_freq - min(DOE_bind_freq)) / (max(DOE_bind_freq) - min(DOE_bind_freq));

% Prepare data for parallel coordinates plot
% Use actual values for vh, np_radius, and koff; normalized values for bind_freq
parallel_data = [reshape(DOE_vh, [], 1), reshape(DOE_np_radius, [], 1), ...
                 reshape(DOE_koff, [], 1), norm_bind_freq];

% Parallel Coordinates Visualization
figure;
parallelcoords(parallel_data, 'labels', {'vh', 'np_radius', 'koff', 'bind_freq'}, ...
               'Color', 0.7*[1 1 1]);
title('Parallel Coordinates Plot of Parameters and Binding Frequency');
grid on;
hold on;

% Highlight high binding frequencies
high_bind_freq_indices = norm_bind_freq > 0.75; % Adjust threshold as needed
subset_parallel_data = parallel_data(high_bind_freq_indices, :);
if ~isempty(subset_parallel_data)
    parallelcoords(subset_parallel_data, ...
                   'labels', {'vh', 'np_radius', 'koff', 'bind_freq'}, ...
                   'Color', 'r', 'LineWidth', 2);
end

hold off;
