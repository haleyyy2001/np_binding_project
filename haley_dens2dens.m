 
%%  
% test Set the parameters
rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;

koff=0.494;

% Generate the surface with clusters
numbPoints = num_clusters;

% Kinetic Params
k0= 0.01;
kon=0.1;

kinetic_params = [k0, kon, koff];

% NP Params
np_rho=0.0002017;
vh = 57;
np_radius = 127;

np_params = [vh, np_radius,np_rho];

figure()


dens_tcr=[];
dens_cluster=[];
prob=[];
i=1;
m=1;

for tcr_per_cluster=1:1:100
    dens_tcr(m)=tcr_per_cluster./ (pi * cluster_radius^2);
    n=1;
    for num_clusters = 1:1:50
 num_tcr=tcr_per_cluster*num_clusters ;

dens_cluster(n)=num_clusters./ (pi * rSurf^2);
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
 [bind_freq, nt_bound, nt_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0)

prob(i)=bind_freq;
i=i+1;

n=n+1;
end
m=m+1;
end
prob=reshape(prob,50,100);
surf(dens_tcr,dens_cluster,prob)
xlabel('density of tcr');ylabel('density of cluster');
zlabel('binding probability')
title([' vh=',num2str(vh),' np radius=',num2str(np_radius),' koff=',num2str(koff) ,' rho=',num2str(np_rho)])
rotate3d



 












%% 
%test
% Set the parameters
rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;

% Generate the surface with clusters
%numbPoints = num_clusters;

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

figure()


dens_tcr=[];
dens_cluster=[];
prob=[];
i=1;
m=1;

for tcr_per_cluster=1:1:100
    dens_tcr(m)=tcr_per_cluster./ (pi * cluster_radius^2);
    n=1;
    for num_clusters = 1:1:50
 num_tcr=tcr_per_cluster*num_clusters ;

dens_cluster(n)=num_clusters./ (pi * rSurf^2);
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
 [bind_freq, nt_bound, nt_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0)

prob(i)=bind_freq;
i=i+1;

n=n+1;
end
m=m+1;
end
prob=reshape(prob,50,100);
surf(dens_tcr,dens_cluster,prob)
xlabel('density of tcr');ylabel('density of cluster');
zlabel('binding probability')
title([' vh=',num2str(vh),' np radius=',num2str(np_radius),' rho=',num2str(np_rho)],'String')
rotate3d




%%
% clusters 


% Set the parameters
rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;

% Generate the surface with clusters
numbPoints = num_clusters;

% Kinetic Params
k0= 0.01;
kon=0.1;
koff=4.94;

kinetic_params = [k0, kon, koff];

% NP Params
np_rho=2.017/10000;
vh = 57;
np_radius = 127;

np_params = [vh, np_radius,np_rho];

figure()


dens_tcr=[];
dens_cluster=[];
prob=[];
i=1;
m=1;

for tcr_per_cluster=1:1:100
    dens_tcr(m)=tcr_per_cluster./ (pi * cluster_radius^2);
    n=1;
    for num_clusters = 1:1:50
 num_tcr=tcr_per_cluster*num_clusters ;

dens_cluster(n)=num_clusters./ (pi * rSurf^2);
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
 [bind_freq, nt_bound, nt_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0)

prob(i)=bind_freq;
i=i+1;

n=n+1;
end
m=m+1;
end
prob=reshape(prob,50,100);
surf(dens_tcr,dens_cluster,prob)
xlabel('density of tcr');ylabel('density of cluster');
zlabel('binding probability')
title([' vh=',num2str(vh),' np radius=',num2str(np_radius),' rho=',num2str(np_rho)],'String')
rotate3d
%% 
%uniform
% Set the parameters
rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;

% Generate the surface with clusters
numbPoints = num_clusters;

% Kinetic Params
k0= 0.01;
kon=0.1;
koff=0.0532;

kinetic_params = [k0, kon, koff];

% NP Params
np_rho=0.0026;
vh =104;
np_radius = 67;

np_params = [vh, np_radius,np_rho];

figure()


dens_tcr=[];
dens_cluster=[];
prob=[];
i=1;
m=1;

for tcr_per_cluster=1:1:100
    dens_tcr(m)=tcr_per_cluster./ (pi * cluster_radius^2);
    n=1;
    for num_clusters = 1:1:50
 num_tcr=tcr_per_cluster*num_clusters ;

dens_cluster(n)=num_clusters./ (pi * rSurf^2);
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
 [bind_freq, nt_bound, nt_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0)

prob(i)=bind_freq;
i=i+1;

n=n+1;
end
m=m+1;
end
prob=reshape(prob,50,100);
surf(dens_tcr,dens_cluster,prob)
xlabel('density of tcr');ylabel('density of cluster');
zlabel('binding probability')
title([' vh=',num2str(vh),' np radius=',num2str(np_radius),' rho=',num2str(np_rho)],'String')
rotate3d












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
function pd = pdf_CovTCR(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius)
    
    tcr_rho = tcr_per_cluster ./ (pi * cluster_radius^2);
    free_rho = (num_tcr - tcr_per_cluster * num_clusters) / (pi* rSurf^2);

    lambda = @(x) ave_CovTCR(x, cluster_radius, np_radius, tcr_rho, free_rho);
    lambda_parents = num_clusters / (pi * (rSurf-cluster_radius)^2);
    Ray_scale = 1 / sqrt(2*lambda_parents*pi);

    pd = [];

    for k = 0:100%number of covered tcr %x is distance 
        fun = @(x, k) (x .* lambda(x) .^ k .* exp(-lambda(x)-x.^2 ./(2 * Ray_scale .^2))) / (Ray_scale.^2 * factorial(k));
        prob = integral(@(x)fun(x,k), 0, 10000);
        pd = [pd, prob];
    end
    xxx=linspace(0,100,101);
   % figure()
    pd = pd ./ sum(pd);
   % plot(xxx,pd)
end

function lambda = ave_CovTCR(x, rnc, rnp, tcr_rho, free_rho)
    A = overlap_area(x, rnc, rnp);
    lambda = A .* tcr_rho + (pi * rnp^2 - A).*free_rho;
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