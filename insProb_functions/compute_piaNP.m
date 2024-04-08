function y = compute_piaNP(r, rctc, ntcr,dtcr)

    % compute the insertion probability of a NP on a nanocluster
    rnc = 100;                  % Originally set to 50
    if nargin<2
        rctc = r/sqrt(2);
        ntcr = 100;             % Originally set to 100
        dtcr = 0.2;              % Originally set to 0.2
    elseif nargin==2
        ntcr = 100;
        dtcr = 0.2;
    end
    
    pia = 1;
    i = 1;
    nsim = 50;
    all = NaN(nsim,ntcr);
    pia_mat = 1;
    while pia >.01 && i<ntcr
        
        for j = 1:nsim
            % Determine initial positions for TCR
            pos_tcr     = rand_scatter_disk(rnc,ntcr,dtcr); 
            % Sample positions for occupied TCR
            pos_np  = sample_min_dist(pos_tcr,i,2*r);           %% Distance between NPs is 2*radius of NP
            if ~isempty(pos_np)
                pos_np  = MC_run(pos_tcr,pos_np,50,2*r,2*r);       %%neighborhood set to 2*r
                expo = exposed_after_binding(pos_tcr,pos_np,2*r);
                all(j,i) = length(expo(1,:))/ntcr;
            end
        end
        pia = mean(all(~isnan(all(:,i)),i));
        pia_mat = [pia_mat,pia];
        all;
        i = i+1;
    end
    x = all(:,~isnan(all(1,:)));
    y = pia_mat(pia_mat>0);

end