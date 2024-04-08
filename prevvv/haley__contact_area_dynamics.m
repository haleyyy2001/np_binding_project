
%% Comprehensive Dynamics of the Contact Area Model
%
close all
%
%
%
%% Contact Area Parameters

nt = 300;                                                                    % Number of TCRs
vh =50;                                                                    % Number of pMHC
N = min(nt,vh);                                                             % Max states of CAM
% TCR Params (Clusters)
vh=50;
np_radius = 100;
% Kinetic Params
k0= 0.01;
kon=0.1;
% NP Params
np_rho=1;
koff=0.1;


Klink = kon/koff;                                                              % pMHC-TCR Dissociation Constant
Ka = k0/koff;





%% Effect of tcr nanocluster radius on bound nanoparticles
% Compare the number of bound NPs for different T cell surface
% organizations. By controlling the radius of TCR nanoclusters using the
% Gaussian filter, we can control the degree of clustering for a fixed number of TCRs.
% (Smaller radius --> More clustered)


plt=true;
cluster_radius = 50;
  rSurf = 1000;
    rTCR = 0.5;
    tcr_num = 300;

   % NP properties
    rNP = 100;
    np_num=50;
    np_rho=1;

%% 

n=20;
    % Surface geometry
    num_clusters = 300/n;

 
 


xx0=0; yy0=0; %centre of disk
 
areaTotal=pi*rSurf^2; %area of disk
  
%Point process parameters
lambda=100; %intensity (ie mean density) of the Poisson process
  
%Simulate Poisson point process
r=rSurf;
numbPoints= num_clusters ;
theta=2*pi*(rand(numbPoints,1)); %angular coordinates
rho=r*sqrt(rand(numbPoints,1)); %radial coordinates
figure()
p = nsidedpoly(1000, 'Center', [xx0, yy0], 'Radius', rSurf);
plot(p,FaceAlpha=0.02)
hold on
%Convert from polar to Cartesian coordinates
[xx,yy]=pol2cart(theta,rho); %x/y coordinates of Poisson points
 
%Shift centre of disk to (xx0,yy0)
xx=xx+xx0;
yy=yy+yy0;
pos=[xx,yy];

pos
t=0;
for i=2:length(pos(:,1))
    if i>length(xx)-t
        break;
    end
    for j=1:i-1
        j
        i
a=dist(pos(j,:),pos(i,:));
if(a(1,1)<=100)
    pos(i,:)=[];
   nn= length(pos(:,1))
    
    i=i-1;
    t=t+1
    continue
end
end


end    

length(pos)


scatter(pos(:,1), pos(:,2));
hold on
for i=1:length(pos(:,1))
%Plotting





p = nsidedpoly(1000, 'Center', [pos(i,1), pos(i,2)], 'Radius', cluster_radius );
plot(p,FaceAlpha=0)
hold

%Simulate Poisson point process
r=cluster_radius;
numbPoints=   n ;
theta=2*pi*(rand(numbPoints,1)); %angular coordinates
rho=r*sqrt(rand(numbPoints,1)); %radial coordinates
 
%Convert from polar to Cartesian coordinates
[xxx,yyy]=pol2cart(theta,rho); %x/y coordinates of Poisson points
 
%Shift centre of disk to (xx0,yy0)
xxx=xxx+pos(i,1);
yyy=yyy+pos(i,2);
posi=[xxx,yyy];
hold

scatter(xxx,yyy);


i=i+1;    
end
xlabel('x');ylabel('y');
axis square;









%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating random positions of points on disk
function pos = generate_pos(rSurf, dmin, num_points);
    rho = rand(1,num_points) + rand(1,num_points);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_points);

    pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];

    i=2;

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
        
         
        i=i+1;
    end
end
%
%
%
%% Generate random clusters of points based on a bimodal distribution
% selection criterion.
function pos = bimodal_pos(rSurf, dmin, num_points);
    init_pts = 100;

    rho = rand(1,100) + rand(1,100);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,100);

    pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];

    i=2;

    while i<length(pos)
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
        i=i+1;
    end
    
    while length(pos) < num_points
        new_tcr = pos(:,randi(length(pos)));
        
        while min(dist(new_tcr', pos)) < 2*dmin
            if rand(1,1) < 0.70
                theta = 2*pi*rand(1,1);
                new_tcr = new_tcr + 2*dmin .* [cos(theta); sin(theta)];
            else
                rho = rand(1,1)+rand(1,1);
                rho(rho>1) = 2-rho(rho>1);
                th = 2*pi*rand(1,1);
                new_tcr = [rSurf*rho*cos(th); rSurf*rho*sin(th)];
            end
        end
        
        pos = [pos, new_tcr];
    end
    
end
%
%
%
%
%% Generate clusters of TCR positions based on gaussian filtering
function pos = generate_clusters(pos, rSurf, num_clusters, cluster_radius, num_tcrs, epsilon)
    % Epsilon controls the base frequency of points not belonging to any
    % cluster. 
    % Choose random centers for Gaussian clusters
    rho = rand(1,num_clusters) + rand(1,num_clusters);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_clusters);

    g_centers = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];
    
    i=1;
    
    while i < length(pos)
        d = min(dist(pos(:,i)',g_centers));
        prob = normpdf(d,0,cluster_radius)./normpdf(0,0,cluster_radius)+epsilon;
        if rand(1,1) > prob
            pos(:,i) = [];
            i=i-1;
        end
        i=i+1;
    end
    
    pos = pos(:,1:num_tcrs);
end
%
%
%
%
%% Compute the number of covered TCRs for each NP
function nt = covered_tcrs(tcr_pos, np_pos, rNP)
    nt = zeros(1,length(np_pos));
    for i = 1:length(np_pos)
        d = dist(np_pos(:,i)',tcr_pos);
        nt(i) = length(d(d<rNP));
    end
end
%
%
%
%% Compute the avidity of a given NP to the T cell surface
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
%
%
%
%
%% Use DBSCAN to perform clustering on TCR positions and yield cluster labels.
function labels = DB_SCAN(pos, minpts)
    kD = pdist2(pos,pos,'euc','Smallest',minpts);
    
    figure()
    plot(sort(kD(end,:)));
    title('k-distance graph')
    xlabel(['Points sorted with ',num2str(minpts),'th nearest distances'])
    ylabel([num2str(minpts),'th nearest distances'])
    grid
    
    epsilon = 40;
    
    labels = dbscan(pos,epsilon,minpts);
    
    figure()
    numGroups = length(unique(labels));
    gscatter(pos(:,1),pos(:,2),labels,hsv(numGroups));
    title(['epsilon = ',num2str(epsilon),' and minpts = ',num2str(minpts)])
    grid
    
    sizeGroups = [];
    
    for  i = 1:numGroups
        sG = length(labels(labels == i));
        sizeGroups = [sizeGroups, sG];
    end
    
    figure()
    histogram(sizeGroups, 10)
    xlabel('TCRs per Nanocluster')
    ylabel('Count')
end

%% Cell Surface Binding
function [K_np, np_pos, nt] = NP_surface_binding(np_params, tcr_params, kinetic_params, plt)
% NP Params = [ NP Radius, Effective Valence, Number of NP Positions, NP concentration]
% TCR Params = [ TCR Radius, TCR Density, Number of Clusters, Radius of Clusters]
% Kinetic Params = [ k0, kon, koff]
    
    rSurf = 1000;
    
    % Kinetic Parameters
    k0 = kinetic_params(1);
    kon = kinetic_params(2);
    koff = kinetic_params(3);
    
    % TCR Properties
    rTCR = tcr_params(1);
    tcr_density = tcr_params(2) / (50^2);
    tcr_num = round(tcr_density * rSurf^2);
    num_clusters = tcr_params(3);
    cluster_radius = tcr_params(4);
    epsilon = tcr_params(5);

    % NP Properties
    rNP = np_params(1);
    vh = np_params(2);
    np_num= np_params(3);
    np_rho= np_params(4);
    
    % Generate positions of TCRs, Clusters and NPs. Estimate covered TCRs.
    tcr_pos = generate_pos(rSurf,rTCR,20*tcr_num);
    tcr_pos = generate_clusters(tcr_pos, rSurf, num_clusters, cluster_radius,tcr_num, epsilon);
    disp('TCR positions generated...')
    np_pos = generate_pos(rSurf, rNP, np_num);
    disp('NP positions generated...')
    nt = covered_tcrs(tcr_pos,np_pos,rNP);
    disp('Covered TCRs computed...')

    % Compute avidity of each NP and probability of binding.
    K_np = [];
    p_bound = [];
    for i=1:size(np_pos,2)
        avidity = np_avidity(vh,nt(i),k0,kon,koff);
        K_np = [K_np, avidity];
        p_bound = [p_bound, np_rho.*avidity /(1 + np_rho.*avidity)];
    end
    disp('Avidity determined...')
    x = rand(1,length(np_pos));
    np_pos = np_pos(:,x<p_bound);
    nt = nt(x<p_bound);
    K_np = K_np(x<p_bound);
    angles = linspace(0,2*pi,500);


    %DB_SCAN(tcr_pos',5)
    %disp('DBSCAN Completed...')

    % Plotting Monte Carlo Sims
    if plt == true
        figure()
        plot(rSurf*cos(angles), rSurf*sin(angles),'k-','Linewidth',2)
        hold on
        for i = 1:length(tcr_pos)
            plot(rTCR*cos(angles)+tcr_pos(1,i),rTCR*sin(angles)+tcr_pos(2,i),'b-')
        end
        for i = 1:size(np_pos,2)
            plot(rNP*cos(angles)+np_pos(1,i), rNP*sin(angles)+np_pos(2,i),'r-')
        end
        
        h = zeros(3, 1);
        h(1) = plot(nan,nan,'or');
        h(2) = plot(nan,nan,'ob');
        h(3) = plot(nan,nan,'ok');
        
        xlim([-rSurf rSurf]);
        ylim([-rSurf rSurf]);
        dim = [0.1 0.6 0.3 0.3];
        str = {['Bound NPs: ',num2str(length(np_pos))],['Average number of TCRs covered per NP: ',num2str(mean(nt(nt~=0)))],['Effective valence of NPs: ',num2str(vh)]};
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
        legend(h,str,'location','southoutside','Fontsize', 12)
        axis off
    end
end