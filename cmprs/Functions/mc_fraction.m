%% 
rSurf = 1000;
cluster_radius = 50;
tcr_radius = 5;
num_tcr=600; 
figure_counter = 0; % Initialize a counter to keep track of the number of figures created
capTCRs=[];
nnn=0;
for num_clusters=30
   for tcr_per_cluster= 0:20
       nnn=nnn+1
 mc_bdtcr=[] ;
 tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
    for ipp=1:500
        mc_bdtcr(ipp)=calculate_binding_capacity(tcr_params)
    end
   capTCRs(nnn)=mean(mc_bdtcr)
   end
end

for np_rho=[0.001]%,0.1,1
for np_radius = [10]  %,20,50
    for  vh = [1]%,10,50, 150
 fig = figure(); % Create a new figure for each vh value
    set(fig, 'Position', [100, 100, 1200, 800]); % Set figur
% Kinetic Params
k0= 0.01;
kon=0.1;
dens_tcr=[];
dens_cluster=[];
prob=[];
prop=[];
 
 pop = 1;
for koff=[0.1]%,0.5,1,5,10,50
i=1;
m=1;
for num_clusters=30
    n=1;
   for tcr_per_cluster= 0:20
 mc_bdtcr=[] ;
  np_params = [vh, np_radius,np_rho];

  tcr_params = [rSurf, num_clusters, cluster_radius, tcr_per_cluster, num_tcr];
kinetic_params = [k0, kon, koff];
[K_np,p_bound, nt, BoundTCRs] = MC_Surface(np_params, kinetic_params, tcr_params, 0);
prob(i)= sum(BoundTCRs)/capTCRs(i)
prop(i)=tcr_per_cluster*num_clusters/600;

i=i+1;
    end
end

subplot(2, 3, pop)
        pop = pop + 1;
        
        plot(prop, prob, '-o')
        
        title(['npradius=', num2str(np_radius), ' koff=', num2str(koff)])
end      
        
        
      
            % Set common x-label and y-label for all subplots
          xlabel('Proportion of tcr in clusters')
ylabel('Boundtcr/Capacity')
            
            % Save the current figure
            cd './mc_res/'
            saveas(fig, ['vh=', num2str(vh), ' rho=', num2str(np_rho), '.fig'])
            close(fig) % Close the current figure
            cd ../
           % figure_counter = 0; % Reset the figure counter
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







