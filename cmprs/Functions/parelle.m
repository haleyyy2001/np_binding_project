% Data extraction from Design of Experiments (assuming you've run your previous code)
 

% Evaluate binding frequency for each experiment 
tcr_params = [1000, 15, 50, 20, 300]; % tcr_params_cluster

 
% Define levels for each parameter
vh_levels = [10, 20, 30,  50,  70,   100];
np_radius_levels = linspace(50, 500, 10);
a = linspace(-1, 2, 10);
koff_levels = 10.^a;

% Factorial design
[DOE_vh, DOE_np_radius, DOE_koff] = ndgrid(vh_levels, np_radius_levels, koff_levels);
 

% Evaluate binding frequency for each experiment
DOE_bind_freq = zeros(num_experiments, 1);
tcr_params = [1000, 15, 50, 20, 300]; % tcr_params_cluster

for i = 1: numel(DOE_vh)
    np_params = [DOE_vh(i), DOE_np_radius(i), 1];
    kinetic_params = [0.01, 0.1, DOE_koff(i)];
   % bind_freq =  
   [bind_freq, nt_bound, nt_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0)
    DOE_bind_freq(i)=sum(bind_freq);
end
% [Your existing code for generating DOE_vh, DOE_np_radius, DOE_koff, and DOE_bind_freq]

% Convert the data to a table for parallelplot
data_table = table(DOE_vh(:), DOE_np_radius(:), DOE_koff(:), DOE_bind_freq(:), ...
                   'VariableNames', {'vh', 'np radius', 'koff', 'bind probability'});

% Create a parallel coordinates plot with data normalization
figure;
p = parallelplot(data_table);
 
% Customize the plot if needed
% For example, to set different scales for each axis or adjust other properties 
