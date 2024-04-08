% Parameters
trial = 100;
np_radius_array = [10, 20, 40];
color_array = ['r', 'g', 'b'];

% Initialize containers for results
capTCRs_all = zeros(length(np_radius_array), 21);
varience_all = zeros(length(np_radius_array), 21);
plotx = zeros(length(np_radius_array), 21*trial); % Change to 2D array

% Main loop for each np_radius value
figure();
hold on;

for np_index = 1:length(np_radius_array)
    np_radius = np_radius_array(np_index);
    color = color_array(np_index);
    
    rSurf = 1000;
    cluster_radius = 50;
    tcr_radius = 5;
    num_tcr = 300; 

    capTCRs = zeros(1, 21);
    varience = zeros(1, 21);
    
    nnn = 0;
    qqq = 0;
    
    for num_clusters = 15
        for tcr_per_cluster = 0:20
            nnn = nnn + 1;
            mc_bdtcr = zeros(1, trial);
            tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
            for ipp = 1:trial
                qqq = qqq + 1;
                mc_bdtcr(ipp) = calculate_binding_capacity(tcr_params, np_radius);
                plotx(np_index, qqq) = mc_bdtcr(ipp); % Save data into 2D array
            end
            capTCRs(nnn) = mean(mc_bdtcr);
            varience(nnn) = var(mc_bdtcr);
        end
    end
    
    capTCRs_all(np_index, :) = capTCRs;
    varience_all(np_index, :) = varience;
    
    yyy = linspace(1, 21*trial, 21*trial);
    plot(yyy, plotx(np_index, :), color);
    for i = 1:21 
        xline(i*trial);
    end
    for n = 0:20
        line([n*trial, trial+n*trial], [capTCRs(n+1), capTCRs(n+1)], 'Color', [0,0,0], 'LineWidth', 2);
    end
end

title('Comparison across np radius values');
legend('np radius:', strsplit(num2str(np_radius_array), ' '));
% Plot histograms
edges = linspace(0, max(plotx(:)), 21); % Update edges based on data range
figure();
handles = zeros(2, length(np_radius_array)); % Store histogram and xline handles

for np_index = 1:length(np_radius_array)
    for mnm = 0:20
        subplot(4, 6, mnm + 1);
        handles(1, np_index) = histogram(plotx(np_index, 1 + mnm*trial : trial + mnm*trial), edges, 'FaceColor', color_array(np_index));
        hold on;
        handles(2, np_index) = xline(capTCRs_all(np_index, mnm + 1), 'Color', color_array(np_index));
        
        titleStr = sprintf('rnp=%d tcr=%d var=%.2f', np_radius_array(np_index), mnm, varience_all(np_index, mnm + 1));
        title(titleStr);
        
        xlims = xlim();
        ylims = ylim();
        text_x = xlims(1) + 0.1 * (xlims(2) - xlims(1));
        text_y = ylims(1) + 0.9 * (ylims(2) - ylims(1));
        text(text_x, text_y, sprintf('Var = %.2f', varience_all(np_index, mnm + 1)));
    end
end

% Add legend after the loop
legendNames = strcat('np radius=', arrayfun(@num2str, np_radius_array, 'UniformOutput', false));
legend(handles(1,:), legendNames); % Display legend using histogram handles for different np_radius values

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
