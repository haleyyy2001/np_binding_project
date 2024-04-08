

%% 

for vh=20
    for np_radius=10
        for koff=[0.1,1,10]
 
np_rho=1;

rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;
 


% Generate the surface with clusters


% Kinetic Params
k0= 0.01;
kon=0.1;

kinetic_params = [k0, kon, koff];


np_params = [vh, np_radius,np_rho];


dens_tcr=[];
dens_cluster=[];
prob=[];
prop=[];



%prop(1)=0;
i=1;
m=1;
num_tcr=600;
%prob(1)=NP_surface_binding(np_params, kinetic_params, [rSurf, 600, cluster_radius, 1, num_tcr], 0);
for tcr_per_cluster=20
    dens_tcr(m)=tcr_per_cluster./ (pi * cluster_radius^2);
    n=1;

   for num_clusters = 0:30
 

dens_cluster(n)=num_clusters./ (pi * rSurf^2);
tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];

[K_np,p_bound, nt, BoundTCRs] = MC_Surface(np_params, kinetic_params, tcr_params,0)
prob(i)=(mean(BoundTCRs(BoundTCRs~=0)));
prop(i)=tcr_per_cluster*num_clusters/600;

i=i+1;
 

% %% 
% for qnum_tcr=[300,500,800,1200,1500,3000]
% qprob=[];
% qi=1;
% qm=1;
% 
% qdens_tcr=[];
% qdens_cluster=[];
% for qtcr_per_cluster=qnum_tcr/50:1:100
%     qdens_tcr(qm)=qtcr_per_cluster./ (pi * cluster_radius^2);
%    
%     qnum_clusters = qnum_tcr/qtcr_per_cluster;
%  
% 
% qdens_cluster(qm)=qnum_clusters./ (pi * rSurf^2);
% qtcr_params = [rSurf, qnum_clusters, cluster_radius, qtcr_per_cluster, qnum_tcr];
%  [qbind_freq, qnt_bound, qnt_pd] = NP_surface_binding(np_params, kinetic_params, qtcr_params, 0)
% 
% qprob(qi)=qbind_freq;
% qi=qi+1;
% 
% 
% qm=qm+1;
% end
% pp=plot3(qdens_tcr,qdens_cluster,qprob)
% pp.LineWidth=3;
% hold on
% text(qdens_tcr(20),qdens_cluster(20),1.2*qprob(20),[num2str(qnum_tcr)],'FontSize',14)
% 
% end
% 
% 





 
    end
end
figure()
i=linspace(0,30,31);

plot(prop,prob,'-o')
xlabel('Proportion of tcr in clusters')
ylabel('Binding Probability')
legend("binding probability")
title([' vh=',num2str(vh),' np radius=',num2str(np_radius),' koff=',num2str(koff) ,' rho=',num2str(np_rho)])

end
    end
end
%% 

function [K_np,p_bound, nt, BoundTCRs] = MC_Surface(np_params, kinetic_params, tcr_params, plt)
    %% Binding of NPs on different cell surfaces
    % Compare the number of bound NPs for different T cell surface
    % organizations.

    % Kinetic
    k0 = kinetic_params(1);
    kon = kinetic_params(2);
    koff= kinetic_params(3);

    % NP
    vh = np_params(1);
    rNP = np_params(2);
    np_rho = np_params(3);
    np_num = 50000;

    % TCR Params
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);
    
    rTCR = 5;

    j=1;
    for num_clusters = tcr_params(2)
        

        disp(['Clusters:', num2str(num_clusters)])

        tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr);  
        disp('TCR positions generated...')
        np_pos = generate_pos(rSurf, 0, np_num);
        disp('NP positions generated...')
        nt = covered_tcrs(tcr_pos,np_pos,rNP);
        disp('Covered TCRs computed...')

        K_np = [];
        p_bound = [];
        BoundTCRs = [];
        for i=1:size(np_pos,2)
            [avidity, bTCR] = np_avidity(vh,nt(i),k0,kon,koff);
            K_np = [K_np, avidity];
            BoundTCRs = [BoundTCRs, bTCR];
            p_bound = [p_bound, avidity /(1 + avidity)];
        end
        disp('Avidity determined...')
        disp(size(p_bound))
        disp(size(np_pos))
        
        r = rand([1,length(p_bound)]);
        bind = find(BoundTCRs > 0);
        np_pos = np_pos(:,bind);
        %K_np = K_np(bind);
        TotalBound = sum(BoundTCRs(bind));
        angles = linspace(0,2*pi,500);
        
        boundNPs = length(np_pos) / length(nt);
        pb = mean(p_bound);

        tot_covTCR = sum(nt(BoundTCRs>0));
        
        [nt_bound, nt_pd, aB, bd_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, plt);
        
        figure(11)
        subplot(3,3,3*j-2)
        histogram(nt(nt~=0),'Normalization','pdf')
        hold on
        plot([1:50],nt_pd ./ sum(nt_pd),'r-')
        xlabel('TCRs Covered Pre')
        
        subplot(3,3,3*j-1)
        histogram(nt(BoundTCRs>0),'Normalization','pdf')
        hold on
        %histogram(nt(bind),'Normalization','pdf')
        plot([1:50], nt_bound ./sum(nt_bound),'r-')
        xlabel('TCRs Covered Post')
        
        disp(find(nt(BoundTCRs>0)==1));
        
        subplot(3,3,3*j)
        histogram(BoundTCRs(BoundTCRs>0),'Normalization','pdf')
        hold on
        plot([1:50], bd_pd ./ sum(bd_pd), 'r-')
        xlabel('Bound')

        %DB_SCAN(tcr_pos',5)
        %disp('DBSCAN Completed...')

        % Plotting Surface
        if plt
            figure(20)
            subplot(1,3,j)
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
            h(4) = plot(nan,nan,'');
            h(5) = plot(nan,nan,'');

            xlim([-rSurf rSurf]);
            ylim([-rSurf rSurf]);
            dim = [0.1 0.6 0.3 0.3];
            str = {['Bound NPs: ',num2str(length(np_pos))],['Average number of TCRs covered per NP: ', ...
            num2str(mean(nt(BoundTCRs>0)))],['Total Covered TCRs: ',num2str(tot_covTCR)],...
            ['Average Bound TCRs: ', num2str(mean(BoundTCRs(BoundTCRs~=0)))], ['Total Bound TCRs: ',num2str(TotalBound)]};
            %annotation('textbox',dim,'String',str,'FitBoxToText','on');
            legend(h,str,'location','southoutside','Fontsize', 12)
            axis off
        end
        j=j+1;
    end
    
    set(gcf, 'Position',[500 500 2000 500]);
end

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
        
        if nsims == 500
            pos(:,i)=[];
            i=i-1;
        end
        i=i+1;
    end
end

%
%
%
%% Generate TCR Nanoclusters (Alternative Gaussian Filtering)
function tcr_pos = gen_clusters(rSurf, cluster_radius, num_clusters, tcr_per_cluster, num_tcr)
    rTCR = 5;
if num_clusters==1
 
rSurf= rSurf-cluster_radius;
dmin =2*cluster_radius;
num_points=num_clusters;
 rho = rand(1,num_points) + rand(1,num_points);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_points);

    nc_pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];
end



if num_clusters>1
    nc_pos = generate_pos(rSurf-cluster_radius, 2*cluster_radius, num_clusters);
end
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





