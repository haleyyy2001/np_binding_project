function y = compute_pia(r,rctc,vh,dtcr,rnbh)
    % This function computes the insertion probability of a TCR on to a NP
    % on the contact area of a NP and a TCRnc. 
    % - r: radius of the NP
    % - v: valency of the NP

    if nargin<4
        dtcr = 0.2;          % Originally set to 0.2
        rnbh = 0.4; %2 * dtcr;    % Region defining 'neighbor pts', maximizes # of neighbours
    end
    
    pia = 1;                % Initialization  of insertion probability
    i = 1;                  
    nsim = 5000;
    all = NaN(nsim,vh+2);
    pia_mat = 1;
    while pia >.01 && i<vh
        
        for j = 1:nsim
            % Determine initial positions for pMHC
            pos_mhc     = rand_scatter_disk(rctc,vh); 
            % Sample positions for occupied pMHC
            pos_tcr_np  = sample_min_dist(pos_mhc,i,dtcr);
            if ~isempty(pos_tcr_np)
                pos_tcr_np  = MC_run(pos_mhc,pos_tcr_np,200,rnbh,dtcr);%%maybe change neighbourhood
                expo = exposed_after_binding(pos_mhc,pos_tcr_np,dtcr); %% For some reason dtcr=5 before??
                all(j,i) = length(expo(1,:))/vh;
            end
        end
        pia = mean(all(~isnan(all(:,i)),i));
        pia_mat = [pia_mat,pia];
        all;
        i = i+1;
    end
    x = all(:,~isnan(all(1,:)));
    y = pia_mat;
end