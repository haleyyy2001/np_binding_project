function [K_np,p_bound, nt, tot_covTCR] = mc_surface_rand( tcr_params)
   

    % TCR Params
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);
    
    rTCR = 5;

    j=1;
    total_covered_tcr=0;
   num_clusters = tcr_params(2) 
   np_pos = generate_pos(rSurf, 50, np_num);


   total_covered_tcr=total_covered_tcr+1;




        for i=1:size(np_pos,2)
            [avidity, bTCR] = np_avidity(vh,nt(i),k0,kon,koff);
            K_np = [K_np, avidity];
            BoundTCRs = [BoundTCRs, bTCR];
           
        end
         p_bound = nt>0
        disp('Avidity determined...')
        disp(size(p_bound))
        disp(size(np_pos))
        
        r = rand([1,length(p_bound)]);
        bind = find(BoundTCRs > 0);
        np_pos = np_pos(:, p_bound);
        %K_np = K_np(bind);
        TotalBound = sum(BoundTCRs(bind));
        angles = linspace(0,2*pi,500);
        
        boundNPs = length(np_pos) / length(nt);
        pb = mean(p_bound);

        tot_covTCR = sum(p_bound)
        
        %DB_SCAN(tcr_pos',5)
        %disp('DBSCAN Completed...')

end
end
%% Generating random positions of points on disk
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
%% Compute the number of covered TCRs for each NP
function nt = covered_tcrs(tcr_pos, np_pos, rNP)
    nt = zeros(1,size(np_pos,2));
    for i = 1:size(np_pos,2)
        d = dist(np_pos(:,i)',tcr_pos);
        nt(i) = length(d(d<rNP));
    end
end

%% Compute the avidity of a given NP to the T cell surface
function [Kav, Bound_tcr] = np_avidity(vh,nt,k0,kon,koff)
    if nt == 0
        Kav = 0;
        Bound_tcr = 0;
        return
    end
    
    N = min(vh,nt);
    i = [1:N-1];
    f0 = k0*vh*nt;
    f = kon.*(vh-i).*(nt-i);
    f = [f0,f];
    b = [1:N] .* koff;
    
    y = cumprod(f./b);
    Kav = sum(y);
    Y0 = 1 / (1 + Kav);

    Bound_pdf = [Y0, y.*Y0];
    Bound_cdf = cumsum(Bound_pdf);
    
    % Random number of bound TCRs estimated from above distribution
    r = rand(1);
    Bound_tcr = find(Bound_cdf >= r,1,'first')-1;
    
end



