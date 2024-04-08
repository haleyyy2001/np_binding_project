function [nt_bound, nt_pd, ave_Bound, bd_pd] = NP_surface_binding(np_params, kinetic_params, tcr_params, plt)
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
    p_bind = Kav ./ (1 + Kav);

    % Probability of covering i TCRs (after binding)
    nt_bound = nt_pd .* p_bind;
    
    % Probability of NP binding (regardless of avidity)
    bind_freq = sum(nt_bound);
    
    % Average Bound TCRs per NP
    [ave_Bound, bd_pd] = boundTCRs(vh,max(nt),k0,kon,koff,nt_pd);
    
    ave_Bound = ave_Bound * bind_freq;

    if plt == 1
        
        figure()
        % Covered TCRs Distribution
        subplot(311)
        plot(nt, nt_pd, 'r-','Displayname','Pre-Bind')
        hold on
        xline(nt*nt_pd'/sum(nt_pd),'r--','DisplayName', [ '\mu = ', num2str(nt*nt_pd'/sum(nt_pd))])
        plot(nt, nt_bound, 'g-','Displayname','Post-Bind')
        xline(nt*nt_bound'/sum(nt_bound),'g--','DisplayName', ['\mu = ',num2str(nt*nt_bound'/sum(nt_bound))])
        xlabel('TCRs covered by NP')
        ylabel('Probability')
        legend()
        title(['P_{bind} =',num2str(bind_freq)])
        
        % Probability of binding versus avidity
        subplot(312)
        plot(nt, p_bind, 'b-','Displayname','Bind Prob.')
        xlabel('TCRs covered by NP')
        ylabel('Probability')
        legend()
        
        % Bound TCR Distribution
        subplot(313)
        plot(nt, bd_pd, 'b-','Displayname','Bound TCRs','Linewidth',1.5)
        hold on
        xline(nt*bd_pd'/sum(bd_pd),'b--','DisplayName', ['\mu = ',num2str(nt*bd_pd'/sum(bd_pd))])
        xlabel('TCRs Bound by NP')
        ylabel('Probability')
        legend()
    end
end

function [Kav, Bound_pdf, Bound_tcr] = np_avidity(vh,nt,k0,kon,koff)
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

    y = cumprod(f./b);
    Kav = sum(y);
    Y0 = 1 / (1 + Kav);
    
    % Distribution Bound TCR
    Bound_pdf = [Y0, y.*Y0];
    Bound_cdf = cumsum(Bound_pdf);
    
    % Random number of bound TCRs estimated from above distribution
    r = rand(1);
    Bound_tcr = find(Bound_cdf >= r,1,'first')-1;
end

function [aveBound, bTCR_pd] = boundTCRs(vh,nt,k0,kon,koff,nt_bound)
    bTCR_pd=[];
    
    for i = 1:nt                  % Bound TCR

        if i > vh
            z = zeros(max(nt)-i+1,1)';
            bTCR_pd = [bTCR_pd, z];
            break
        end
        p_prod = 0;
        for j = i:nt            % Covered TCR
            p_cov = nt_bound(j) / sum(nt_bound);
            [Kav, Bound_pdf, Bound_tcr] = np_avidity(vh,j,k0,kon,koff);
            p_condi = Bound_pdf(i+1);
            p_prod = p_prod + p_condi * p_cov;
        end

        bTCR_pd = [bTCR_pd, p_prod];
    end

    aveBound = [1:length(bTCR_pd)] * bTCR_pd';

end