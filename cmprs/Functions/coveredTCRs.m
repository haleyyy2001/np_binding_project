function nt = coveredTCRs(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius, plt)
% rSurf = 
% free_rho =  Free tcr density (not in cluster)
    nsamples = 100000;

    %nt_per_nc = poissrnd(tcr_per_cluster,1000,1);
    tcr_rho = tcr_per_cluster ./ (pi * cluster_radius^2);
    free_rho = (num_tcr - tcr_per_cluster * num_clusters) / (pi* rSurf^2);
    rNP = np_radius;

    % Poisson point process
    if num_clusters == 0
        tcr_rho = free_rho;
        num_clusters = 1;   % To avoid error in next line. Does not affect results.
    end
        
    lambda_parents = num_clusters / (pi * (rSurf-cluster_radius)^2);

    % Spherical contact distribution
    x = [0:1000];

    Ray_scale = 1 / sqrt(2*lambda_parents*pi);
    y = raylrnd(Ray_scale,nsamples,1);              % Generate Rayleigh random variable for distance to nearest cluster center
    
    A = overlap_area(y, cluster_radius, rNP);
    lambda = A .* tcr_rho + (pi * rNP^2 - A).*free_rho;
    
    nt = poissrnd(lambda,nsamples,1);
    
    A2 = overlap_area(x, cluster_radius, rNP);
    lambda2 = A2 .* tcr_rho + (pi * rNP^2 - A2).*free_rho;

    if plt
        angles = linspace(0,2*pi,500);

        figure()
        plot(x, lambda2,'r-')
        xlim([0, 2*cluster_radius])
        xlabel('NP Position from cluster center')
        ylabel('Number of covered TCR')

        figure()
        histogram(nt,(0.5:round(max(nt))+0.5),'Normalization','pdf')
        xlabel('Covered TCRs')
        ylabel('Frequency')
    end
end

function area = overlap_area(x, rnc, rnp)
% x =  NP distance from center of NC
% rnc = Radius of nanocluster
% rnp = Radius of nanoparticle

    if rnp > rnc
        tmp = rnc;
        rnc = rnp;
        rnp = tmp;
    end

    h = 1./(2.*x).*sqrt(2.*x.^2.*rnp^2+2.*x.^2.*rnc^2+2*rnp^2*rnc^2-rnp^4-rnc^4-x.^4);
    
    a1 = rnp^2.*asin(h./rnp)+rnc^2.*asin(h./rnc)-x.*h;
    
    a2 = rnp^2 .*asin(h./rnp)-h.*sqrt(rnp^2-h.^2);
    a2 = a2 - rnc^2 .*asin(h./rnc)+ h.*sqrt(rnc^2-h.^2);
    a2 = pi * rnp^2 - a2;
    
    area = a1;
    
    cond = sqrt(rnc^2 -h.^2);
    
    area(x <= cond) = a2(x <= cond);
    
    area = real(area);

end