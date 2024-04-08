for vh=10:30:200
    for np_radius=[10,50,100,250]
        for koff=[0.02,0.05,0.1,0.5,1,5,20,60,100]
close all
np_rho=2.017/10000;




rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;



% Generate the surface with clusters


% Kinetic Params
k0= 0.01;
kon=0.1;

kinetic_params = [k0, kon, koff];


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
hold on




%% 
for qnum_tcr=[300,500,800,1200,1500,3000]
qprob=[];
qi=1;
qm=1;

qdens_tcr=[];
qdens_cluster=[];
for qtcr_per_cluster=qnum_tcr/50:1:100
    qdens_tcr(qm)=qtcr_per_cluster./ (pi * cluster_radius^2);
   
    qnum_clusters = qnum_tcr/qtcr_per_cluster;
 

qdens_cluster(qm)=qnum_clusters./ (pi * rSurf^2);
qtcr_params = [rSurf, qnum_clusters, cluster_radius, qtcr_per_cluster, qnum_tcr];
 [qbind_freq, qnt_bound, qnt_pd] = NP_surface_binding(np_params, kinetic_params, qtcr_params, 0)

qprob(qi)=qbind_freq;
qi=qi+1;


qm=qm+1;
end
pp=plot3(qdens_tcr,qdens_cluster,qprob)
pp.LineWidth=3;
hold on
text(qdens_tcr(20),qdens_cluster(20),1.2*qprob(20),[num2str(qnum_tcr)],'FontSize',14)

end






%% 
xlabel('density of tcr');ylabel('density of cluster');
zlabel('binding probability')
 prob=reshape(prob,1,5000);
 ii = min(find(prob== max(prob)));
pmhc_c=vh/(0.5*pi*(np_radius)^2);
tcr_per_cluster=linspace(1,100,100)

tcr_per_cluster=  tcr_per_cluster(ii/50)
title([' vh=',num2str(vh),' np radius=',num2str(np_radius),' koff=',num2str(koff) ,' rho=',num2str(np_rho),'peak=',num2str(tcr_per_cluster),"pmhc concentration=",pmhc_c])
text(3,5,num2str(pmhc_c))
rotate3d



cd 'HALEY/ac_list'
savefig(['vh=',num2str(vh),' np radius=',num2str(np_radius),' koff=',num2str(koff) ,' rho=',num2str(np_rho),'peak=',num2str(tcr_per_cluster),'.fig'])
cd ../../
    end
end
end


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


