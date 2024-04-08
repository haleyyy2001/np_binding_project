% TCR Params (Clusters)
rSurf = 1000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 20*15;

tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

% Kinetic Params
k0= 0.01;
kon=0.1;
koff=0.1;

kinetic_params = [k0, kon, koff];

% NP Params
np_rho=1;
vh = 50;
np_radius = 100;

np_params = [vh, np_radius,np_rho];

NP_surface_binding(np_params, kinetic_params, tcr_params, 1)

% TCR Params (Uniform Surface)
rSurf = 1000;
num_clusters = 20*15;
cluster_radius = 50;
tcr_per_cluster = 1;
num_tcr = 20*15;

tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

NP_surface_binding(np_params, kinetic_params, tcr_params, 1)



    %% coveredTCRs(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius, plt)
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
  
         poss=[];
        Kav = [];
    for i =[1:round(max(nt))+1]
     
        % Probability of covering i TCRs (before binding):
    nt_pd = pdf_CovTCR(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius);
    nt_pd(1) = [];                                                          % Remove first element corresponding to nt=0
    nt_pd(length(i)+1:end)=[];                                             % Remove elements outside of range
    
        K = np_avidity(vh,i,k0,kon,koff);
        % Probability of binding (as function of avidity)
         p_bind = np_rho*K ./ (1 + np_rho * K);

    % Probability of covering i TCRs (after binding)
     nt_bound = nt_pd .* p_bind;
        prp = nt(i)*nt_bound;
        poss = [poss,prp];
    end
    poss
    
  



%% 
function [bind_freq, nt_bound, nt_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, plt)
    % TCR Params
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);

    % Kinetic Params
    k0= kinetic_params(1);
    kon=kinetic_params(2);
    koff=kinetic_params(3);

    % NP Params
    vh = np_params(1);
    np_radius = np_params(2);
    np_rho = np_params(3);

    nt = [1:50];                                                            % Define range of covered TCRs to consider
    
    % Probability of covering i TCRs (before binding):
    nt_pd = pdf_CovTCR(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius);
    nt_pd(1) = [];                                                          % Remove first element corresponding to nt=0
    nt_pd(length(nt)+1:end)=[];                                             % Remove elements outside of range
    
    % Compute NP avidity
    Kav = [];
    for i = nt
        K = np_avidity(vh,i,k0,kon,koff);
        Kav = [Kav,K];
    end
    
    % Probability of binding (as function of avidity)
    p_bind = np_rho*Kav ./ (1 + np_rho * Kav);

    % Probability of covering i TCRs (after binding)
    nt_bound = nt_pd .* p_bind;
    
    % Probability of NP binding (regardless of avidity)
    bind_freq = sum(nt_bound);

    if plt
        figure()
        subplot(211)
        plot(nt, nt_pd, 'r-','Displayname','Pre-Bind')
        hold on
        plot(nt, nt_bound, 'g-','Displayname','Post-Bind')
        xlabel('TCRs covered by NP')
        ylabel('Probability')
        legend()
        title(['Cluster Landscape: P_{bind} =',num2str(bind_freq)])

        subplot(212)
        plot(nt, p_bind, 'b-','Displayname','Bind Prob.')
        xlabel('TCRs covered by NP')
        ylabel('Probability')
        legend()
    end
end

function Kav = np_avidity(vh,nt,k0,kon,koff)
    if nt == 0
        Kav = 0;
        return
    end
    
    N = min(vh,nt);
    i = [1:N-1];
    f0 = k0*vh*nt;
    f = kon.*(vh-i).*(nt-i);
    f = [f0,f];
    b = [1:N] .* koff;
Kav = cumprod(f./b);
    Kav = sum(Kav);
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














    