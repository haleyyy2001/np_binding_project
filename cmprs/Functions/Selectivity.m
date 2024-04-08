function [selectivity,p_bound] = Selectivity(np_params,kinetic_params, tcr_params, plt)

    % Kinetic Params
    k0= kinetic_params(1);
    kon= kinetic_params(2);
    koff= kinetic_params(3);

    % NP Params
    vh = np_params(1);
    np_radius = np_params(2);
    np_rho= np_params(3);

    % TCR Params
    rSurf = tcr_params(1);
    num_clusters = tcr_params(2);
    cluster_radius = tcr_params(3);
    tcr_per_cluster = tcr_params(4);
    num_tcr = tcr_params(5);

    nt = [1:100];

    K_np = [];
    p_bound = [];
    selectivity = [];

    for i=nt
        avidity = np_avidity(vh,i,k0,kon,koff);
        K_np = [K_np, avidity];
        p_bound = [p_bound, avidity /(1 + avidity)];

        if i>1
            alpha = log10(p_bound(i))-log10(p_bound(i-1));
            alpha = alpha / (log10(i)-log10(i-1));
            selectivity = [selectivity,alpha];
        end
    end
    
    x = nt / (pi * np_radius^2);
    
    x1 = tcr_per_cluster / (pi * cluster_radius^2);
    x2 = num_tcr / (pi * rSurf^2);
    
    if plt == 1
        figure()
        subplot(211)
        plot(x, p_bound)
        xline(x1,'g--','DisplayName',['Cluster Density = ',num2str(x1)])
        xline(x2,'r--','DisplayName',['Average Density = ',num2str(x2)])
        xlabel('TCR Density')
        ylabel('Surface Coverage')
        
        subplot(212)
        plot(x(1:end-1), selectivity)
        xline(x1,'g--','DisplayName',['Cluster Density = ',num2str(x1)])
        xline(x2,'r--','DisplayName',['Average Density = ',num2str(x2)])
        xlabel('TCR Density')
        ylabel('Selectivity')
        
        legend()
    end
end
%--------------------------------------------------------------------------
%% Compute the avidity of a given NP to the T cell surface
%--------------------------------------------------------------------------
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