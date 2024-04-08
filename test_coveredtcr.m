
% Kinetics
kon = 0.1;
k0 = 0.01;
koff= 0.1;

% TCRs
rTCR = 5;
tcd = 1;
num_clusters= 15;
cluster_radius = 50;

% NPs
rNP = 100;
vh = 1;
np_num = 50;
tcr_per_cluster=2;


rSurf = 1000;
num_clusters = 15;
cluster_radius = 50;
tcr_per_cluster = 20;
num_tcr = 20*15;

np_radius=rNP;
plt=true;
num_tcr=num_clusters;
coveredTCRs(rSurf, num_tcr, num_clusters, cluster_radius, tcr_per_cluster, np_radius, plt)