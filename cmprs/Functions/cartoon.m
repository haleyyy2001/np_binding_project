% Initialization
trial = 1000;
np_radius = 50;
rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;
num_tcr = 600;
capTCRs = [];
nnn = 0;
qqq = 0;

% Specify figure dimensions
fig = figure('Position', [0, 0, 1120, 840]);

writerObj = VideoWriter('distribution_shift_video.avi');
writerObj.FrameRate = 30;  % Set FrameRate
open(writerObj);

%% Main Loop
for num_clusters = 30
    for tcr_per_cluster = 0:20
        nnn = nnn + 1;
        mc_bdtcr = [];
        tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
        
        for ipp = 1:trial
            qqq = qqq + 1;
            mc_bdtcr(ipp) = calculate_binding_capacity(tcr_params, np_radius);
            plotx(qqq) = mc_bdtcr(ipp);
        end
        
        capTCRs(nnn) = mean(mc_bdtcr);
        
        % Clear previous plots
        clf;
        
        % Create Histogram
        edges = 0:1:max(plotx);
        h = histogram(plotx(1 + (nnn - 1) * trial : nnn * trial), edges);
        
        % Get peak value
        [~, peakIdx] = max(h.Values);
        peakValue = h.BinEdges(peakIdx) + h.BinWidth / 2;
        
        % Add vertical line marking the mean
        xline(capTCRs(nnn), 'Color', 'r', 'Label', ['Mean = ' num2str(capTCRs(nnn))]);
        
        % Add vertical line marking the peak
        xline(peakValue, 'Color', 'b', 'Label', ['Peak = ' num2str(peakValue)]);
        
        % Add title
        title(['TCRs per cluster = ', num2str(tcr_per_cluster)]);
        
        % Capture the current figure and write to video
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end
end

% Close video writer
close(writerObj);

%% Other plotting logic and functions remain unchanged...
% ...









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

%



