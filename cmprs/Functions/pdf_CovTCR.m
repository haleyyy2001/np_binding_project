function pd = pdf_CovTCR(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius)
    
    % Densities of TCRs within (tcr_rho) and between (free_rho)
    % nanoclusters
    tcr_rho = tcr_per_cluster ./ (pi * cluster_radius^2);
    free_rho = (num_tcr - tcr_per_cluster * num_clusters) / (pi* rSurf^2);
    
    if num_clusters == 0
        lambda = free_rho * pi * np_radius^2;
        pd = poisspdf([0:100], lambda);
        pd = pd ./ sum(pd);
        return
    end
    
    % Poisson mean of TCRs covered by NP landing a distance x from
    % Nanocluster.
    lambda = @(x) ave_CovTCR(x, cluster_radius, np_radius, tcr_rho, free_rho);
    
    % Poisson mean density of nanoclusters on cell surface.
    lambda_parents = num_clusters / (pi * (rSurf-cluster_radius)^2);
    
    % Rayleigh scale parameter for expected distance to nearest
    % nanocluster.
    Ray_scale = 1 / sqrt(2*lambda_parents*pi);

    pd = [];

    % Integral of Poisson and Rayleigh distributions to calculate Poisson
    % pdf of TCRs covered by NP.
    for k = 0:100
        fun = @(x, k) (x .* lambda(x) .^ k .* exp(-lambda(x)-x.^2 ./(2 * Ray_scale .^2))) / (Ray_scale.^2 * factorial(k));
        prob = integral(@(x)fun(x,k), 0, 10000);
        pd = [pd, prob];
    end
    
    pd = pd ./ sum(pd);
end

%% Local Functions
% Expected TCRs covered by NP a distance x from nearest NC
function lambda = ave_CovTCR(x, rnc, rnp, tcr_rho, free_rho)
    A = overlap_area(x, rnc, rnp);
    lambda = A .* tcr_rho + (pi * rnp^2 - A).*free_rho;
end

% Quantify area of overlap between NP and NC
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