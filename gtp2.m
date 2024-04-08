
% Set the parameters
rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;
tcr_per_cluster = 30;
num_clusters = round(300/tcr_per_cluster);

xx0 = 0; yy0 = 0; % Centre of disk

% Generate the surface with clusters
numbPoints = num_clusters;
theta = 2*pi*(rand(numbPoints,1)); % Angular coordinates
rho = rSurf*sqrt(rand(numbPoints,1)); % Radial coordinates
[xx, yy] = pol2cart(theta,rho); % x/y coordinates of Poisson points
xx = xx + xx0; yy = yy + yy0; % Shift centre of disk to (xx0,yy0)
cluster_centers = [xx, yy];

% Check for overlaps in clusters and remove
overlap = true;
while overlap == true
    overlap = false;
    for i = 1:length(cluster_centers)
        for j = (i+1):length(cluster_centers)
            if norm(cluster_centers(i,:) - cluster_centers(j,:)) <= 2*cluster_radius
                cluster_centers(j,:) = [];
                overlap = true;
                break;
            end
        end
        if overlap == true
            break;
        end
    end
end

figure
p = nsidedpoly(1000, 'Center', [xx0, yy0], 'Radius', rSurf);
plot(p,FaceAlpha=0.02)
hold on

for i = 1:length(cluster_centers)
    % Plot each cluster
    p = nsidedpoly(1000, 'Center', cluster_centers(i,:), 'Radius', cluster_radius);
    plot(p, 'FaceAlpha', 0)

    % Generate TCRs within each cluster
    numbPoints = tcr_per_cluster;
    theta = 2*pi*(rand(numbPoints,1)); % Angular coordinates
    rho = cluster_radius*sqrt(rand(numbPoints,1)); % Radial coordinates
    [xx, yy] = pol2cart(theta,rho); % x/y coordinates of Poisson points
    xx = xx + cluster_centers(i,1); yy = yy + cluster_centers(i,2); % Shift centre of disk to cluster center
    tcrs = [xx, yy];

    % Plot TCRs
    scatter(tcrs(:,1), tcrs(:,2));
end

xlabel('x'); ylabel('y');
title('tcr per cluster',tcr_per_cluster)
axis square;
