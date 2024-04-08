% TCR Params (Clusters)
close all;
rSurf = 1000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 20*15;

tcr_params_cluster = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

% Kinetic Params
k0= 0.01;
kon=0.1;
% NP Params
np_rho=1;

% TCR Params (Uniform Surface)
rSurf = 1000;
num_clusters = 300;
cluster_radius = 50;
tcr_per_cluster = 1;
num_tcr = 20*15;
tcr_params_uniform = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
vh_max=100;
np_radius_max=500;
int=1;
G1=linspace(1,vh_max,vh_max);
G2=linspace(50,np_radius_max, 20);
a=linspace(-1,2,30);
G3=10.^a;
f_max=-1e+10;
f_min=10000;
nR=5;

%for r=1:nR
% ... [the initialization of the variables remains unchanged] ...














% Pre-allocate arrays for efficiency
DIFC = zeros(100, length(a));
DIFD = zeros(100,length(a));

% Populate the DIFC and DIFD arrays
for i=1:100
    for j=1:length(a)
        DIFC(i, j) = NP_surface_binding([G1(i), 20, np_rho], [k0, kon, G3(j)], tcr_params_cluster, 0);
        DIFD(i, j) = NP_surface_binding([G1(i), 20, np_rho], [k0, kon, G3(j)], tcr_params_uniform, 0);
    end
end

% Plotting
figure()
surf(G1, G3, DIFC')
xlabel('valence'); ylabel('koff');
title('cluster valence v.s. koff; np radius=20 rho=1'); % Adjust title appropriately

figure()
surf(G1, G3, DIFD')
xlabel('valence'); ylabel('koff');
title('uniform valence v.s. koff; np radius=20 rho=1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DIFC = zeros(20,length(a));
DIFD = zeros(20, length(a));

% Populate the DIFC and DIFD arrays
for i=1:20
    for j=1:length(a)
          DIFC(i, j) = NP_surface_binding([5, G2(i), np_rho], [k0, kon, G3(j)], tcr_params_cluster, 0);
        DIFD(i, j) = NP_surface_binding([5, G2(i), np_rho], [k0, kon, G3(j)], tcr_params_uniform, 0);
    end
end

% Plotting
figure()
surf(G2, G3, DIFC')
xlabel('np radius'); ylabel('koff');
title('cluster np radius v.s. koff; valence=5 rho=1');
figure()
surf(G2, G3, DIFD')
xlabel('np radius'); ylabel('koff');
title('uniform np radius v.s. koff; valence=5 rho=1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

DIFC = zeros(100, 20);
DIFD = zeros(100, 20);


% Populate the DIFC and DIFD arrays
for i=1:100
    for j=1:20
        DIFC(i, j) = NP_surface_binding([G1(i), G2(j), np_rho], [k0, kon, 0.1], tcr_params_cluster, 0);
        DIFD(i, j) = NP_surface_binding([G1(i), G2(j), np_rho], [k0, kon, 0.1], tcr_params_uniform, 0);
    end
end

% Plotting
figure()
surf(G1, G2, DIFC')
xlabel('valence'); ylabel('np radius');
title('cluster valence v.s. np radius; koff=0.1 rho=1');


figure()
surf(G1, G2, DIFD')
xlabel('valence'); ylabel('np radius');
title('uniform valence v.s. np radius; koff=0.1 rho=1');


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
