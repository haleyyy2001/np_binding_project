% Design of Experiments for Parameter Sensitivity

% Define levels for each parameter
vh_levels = [1, 50, 100];
np_radius_levels = [50, 275, 500];
koff_levels = [0.01, 5, 10];

% Factorial design
[DOE_vh, DOE_np_radius, DOE_koff] = ndgrid(vh_levels, np_radius_levels, koff_levels);
num_experiments = numel(DOE_vh);

% Evaluate binding frequency for each experiment
DOE_bind_freq = zeros(num_experiments, 1);

tcr_params = [1000, 15, 50, 20, 300]; % tcr_params_cluster

for i = 1:num_experiments
    np_params = [DOE_vh(i), DOE_np_radius(i), 1];
    kinetic_params = [0.01, 0.1, DOE_koff(i)];
    
    bind_freq = NP_surface_binding(np_params, kinetic_params, tcr_params, 0);
    DOE_bind_freq(i) = mean(bind_freq(:));
end

% Fit a linear model using formula
formula = 'DOE_bind_freq ~ DOE_vh + DOE_np_radius + DOE_koff';
lm = fitlm(table(DOE_vh(:), DOE_np_radius(:), DOE_koff(:), DOE_bind_freq(:), ...
        'VariableNames', {'DOE_vh', 'DOE_np_radius', 'DOE_koff', 'DOE_bind_freq'}), formula);

% Display the coefficients and p-values
disp(lm.Coefficients)

% Display p-values for each variable
disp('P-values for vh, np_radius, koff:')
disp(lm.Coefficients.pValue(2:end))

