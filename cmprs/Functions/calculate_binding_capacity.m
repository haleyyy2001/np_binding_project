

function [total_binded_np ] = calculate_binding_capacity(tcr_params)

    % Given parameters
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);

    rTCR = 5;
    nanoparticle_radius = 50;

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
