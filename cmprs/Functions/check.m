% Parameters
trial = 1000 ;
np_radius_array = [10, 20, 40];  % Different np_radius values
color_array = [[0.6, 0.8, 1]; [1, 0.7, 0.9]; [1, 0.8, 0.6]];  % Light blue, light pink, light orange

% Initialize containers for results
plotx = zeros(length(np_radius_array), 21 * trial); 

% Main loop for each np_radius value
figure();
hold on;

plot_handles = zeros(1, length(np_radius_array));  % Array for storing plot handles

for np_index = 1:length(np_radius_array)
    np_radius = np_radius_array(np_index);
    color = color_array(np_index, :);
    capTCRs = zeros(1, 21);
    varience = zeros(1, 21);
    qqq = 0;
    
    for tcr_per_cluster = 0:20
        mc_bdtcr = zeros(1, trial);
        tcr_params = [1000, 15, 50, tcr_per_cluster, 300];
        for ipp = 1:trial
            qqq = qqq + 1;
            mc_bdtcr(ipp) = calculate_binding_capacity(tcr_params, np_radius);
            plotx(np_index, qqq) = mc_bdtcr(ipp);
        end
        capTCRs(tcr_per_cluster + 1) = mean(mc_bdtcr);
        varience(tcr_per_cluster + 1) = var(mc_bdtcr);
    end
    
    yyy = linspace(1, 21 * trial, 21 * trial);
    plot_handles(np_index) = plot(yyy, plotx(np_index, :), 'Color', color, 'LineWidth', 1.5);  % Storing the handle
    for i = 1:21
        xline(i * trial, '--k', 'LineWidth', 0.5);
    end
    for n = 0:20
        line([n*trial, trial+n*trial], [capTCRs(n+1), capTCRs(n+1)], 'Color', [0, 0, 0]);
    end
end

xlabel('Trials');
ylabel('Binding Capacity');
title('Comparison across np radius values');
legend(plot_handles, arrayfun(@(x) ['np radius = ', num2str(x)], np_radius_array, 'UniformOutput', false));

set(gca, 'FontSize', 14); % Setting font size for readability



% Histogram plots
for np_index = 1:length(np_radius_array)
    figure();
    for mnm = 0:20
        subplot(4, 6, mnm + 1);
        binEdges = 0:5:max(plotx(np_index, :)) + 5;  % example bin edges from 0 to max+5 with 10 units apart
histogram(plotx(np_index, 1 + mnm * trial : trial + mnm * trial), 'BinEdges', binEdges, 'FaceColor', color_array(np_index, :));

        hold on;
        xline(mean(plotx(np_index, 1 + mnm * trial : trial + mnm * trial)), 'Color', color_array(np_index, :), 'LineWidth', 1.5);
        
        current_mean = mean(plotx(np_index, 1 + mnm * trial : trial + mnm * trial));
        current_var = var(plotx(np_index, 1 + mnm * trial : trial + mnm * trial));
        textStr = sprintf('Mean=%.2f\nVar=%.2f', current_mean, current_var);
        xlims = xlim();
        ylims = ylim();
        text_x = 5;
        text_y = 300;
        text(text_x, text_y, textStr, 'FontSize', 8);
        
        xlim([0, max(plotx(np_index, :)) *1.2]);
        ylim([0, 400]);
        
        set(gca, 'FontSize', 10);
    end
    sgtitle(['np radius = ', num2str(np_radius_array(np_index))], 'FontSize', 14);
end

% Your function implementations remain unchanged...

% Your function implementations remain unchanged
function [total_binded_np ] = calculate_binding_capacity(tcr_params,np_radius)

    % Given parameters
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);

    rTCR = 5;
    nanoparticle_radius = np_radius;

    % Generate initial TCR positions
    tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr); % Assuming generate_pos returns a 2xN array.

    % Initialize variables
    total_binded_np = 0;
    nanoparticle_positions = [];

    for i = 1:num_tcr
        current_tcr = tcr_pos(:, i);

        % Check if nanoparticle can be placed without overlap
        is_overlap = false;
        for j = 1:size(nanoparticle_positions, 2)
            existing_np = nanoparticle_positions(:, j);
            distance = norm(current_tcr - existing_np);
            
            if distance < 2 * nanoparticle_radius
                is_overlap = true;
                break;
            end
        end

        % Place nanoparticle if no overlap
        if ~is_overlap
            nanoparticle_positions = [nanoparticle_positions, current_tcr];
            total_binded_np = total_binded_np + 1;
        end
    end

    % Display the total number of binded nanoparticles
    disp(['Total binded nanoparticles: ', num2str(total_binded_np)]);
end

function pos = generate_pos(rSurf, dmin, num_points);
    rho = rand(1,num_points) + rand(1,num_points);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_points);

    pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];

    i=3;
  
if length(pos)<=2 
   
    return 
else
    while i<=length(pos)
        d = dist(pos(:,i)',pos(:,1:i-1));
        min_d = min(d);
        nsims = 0;
        while min_d < 2*dmin && nsims < 500
            rho = rand(1,1)+rand(1,1);
            rho(rho>1) = 2-rho(rho>1);
            th = 2*pi*rand(1,1);
            pos(:,i) = [rSurf*rho*cos(th); rSurf*rho*sin(th)];
        
            d = dist(pos(:,i)',pos(:,1:i-1));
            min_d = min(d);
            nsims = nsims+1;
        end
        
        if nsims == 500
            pos(:,i)=[];
            i=i-1;
        end
        i=i+1;
    end
end
end

%
%
%
%% Generate TCR Nanoclusters (Alternative Gaussian Filtering)
function tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr)
    rTCR = 5;
    nc_pos = generate_pos(rSurf-cluster_radius, 2*cluster_radius, num_clusters);
    tcr_pos = [];




    for i = 1:num_clusters
        lambda = poissrnd(tcr_per_cluster);
        %lambda = tcr_per_cluster;
        temp_pos = generate_pos(cluster_radius, rTCR, lambda);
        temp_pos = temp_pos + nc_pos(:,i);
        tcr_pos = [tcr_pos, temp_pos];
    end
    
    free_tcr = num_tcr -length(tcr_pos);
    
    if free_tcr < 0
        tcr_pos = tcr_pos(:,1:num_tcr);
    else
        temp_pos = generate_pos(rSurf, rTCR, free_tcr);
        tcr_pos = [tcr_pos, temp_pos];
    end
end

