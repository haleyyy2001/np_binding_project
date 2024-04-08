clear all
%path(path, '/home/louis/Desktop/NP-Revival/Functions/')

%C = 1; U = 2;
choose_np = 1;
choose_surface = 1;

% Kinetic and NP Params

if choose_np == 1  
    k0= 0.01;
    kon=0.1;
    koff=53;
    
    np_rho=0.54;
    vh = 101;
    np_radius = 72;
    
elseif choose_np == 2
    k0= 0.01;
    kon=0.1;
    koff=1.22;
    
    np_rho=1;
    vh = 34;
    np_radius = 115;
    
else
    k0= 0.01;
    kon=0.1;
    koff=4.94;
    
    np_rho=0.02017;
    vh = 28.5;
    np_radius = 127;
end

% TCR Params
if choose_surface == 1
    rSurf = 3000;
    num_clusters = 15;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 20*15;
    tcr_radius = 5;
    
elseif choose_surface == 2
    rSurf = 3000;
    num_clusters = 0;
    cluster_radius = 50;
    num_tcr = 20*15;
    tcr_per_cluster = 20; %num_tcr / (rSurf^2 / cluster_radius^2 - num_clusters);
    tcr_radius = 5;
    
else
    rSurf = 1000;
    num_clusters = 15;
    cluster_radius = 50;
    tcr_per_cluster = 20;
    num_tcr = 20*15;
    tcr_radius = 5;
end

k0 = k0 * np_rho;

kinetic_params = [k0, kon, koff];
np_params = [vh, np_radius, np_rho];
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

% Distribution of covered TCRs and bound TCRs
%[nt_bound, nt_pd, bd_pd, Kav, p_bind] = NP_surface_binding(np_params, kinetic_params, tcr_params, 0);

% Surface Coverage and Selectivity
Selectivity(np_params, kinetic_params, tcr_params, 1);

% MC Surface Coverage
[K_np, p_bound, nt, BoundTCRs] = MC_Surface(np_params, kinetic_params, tcr_params, 1);

figure()
subplot(221)
scatter(K_np, p_bound)
xlabel('Avidity')
ylabel('Bind Prob')

subplot(222)
scatter(nt, BoundTCRs)
xlabel('Covered TCRs')
ylabel('Bound TCRs')

subplot(223)
histogram(BoundTCRs(BoundTCRs~=0))
xlabel('Bound TCRs')

subplot(224)
histogram(nt(nt~=0))
xlabel('Covered TCRs')

%%


cluster_params = [3000, 15, 50, 20, 15*20];
uni_params = [3000,0,50,20,15*20];
kinetic_params = [0.01,0.1,20];

cluster_select = [];
cluster_Bound = [];
uni_select = [];
uni_Bound = [];

x = -4:0.5:0;
x = 10.^x;

for rho = x
    
    np_params = [50, 100, 1];
    kinetic_params = [rho*0.01,0.1,20];
    
    [selectivity, p_bound] = Selectivity(np_params, kinetic_params, cluster_params, 1);
    cluster_select = [cluster_select, find(selectivity == max(selectivity))];
    
    [nt_bound, nt_pd, aveBound, bd_pd] = NP_surface_binding(np_params, kinetic_params, cluster_params, 0);
    cluster_Bound = [cluster_Bound, aveBound];
    
    [nt_bound, nt_pd, aveBound, bd_pd] = NP_surface_binding(np_params, kinetic_params, uni_params, 0);
    uni_Bound = [uni_Bound, aveBound];
end

figure()
plot(x, cluster_select)
xlabel('Valence')
ylabel('Number of TCRs of peak selectivity')

figure()
semilogy(x, cluster_Bound, 'DisplayName','Cluster')
hold on
semilogy(x, uni_Bound, 'DisplayName', 'Uniform')
xlabel('Valence')
ylabel('Bound TCRs')
legend()
