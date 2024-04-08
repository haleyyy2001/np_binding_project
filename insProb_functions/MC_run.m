function y = MC_run(pos,sam,nsim,rngb,rocc)
    % This function allows for a MC loop in the manner exposed by Hlavacek
    % et. al.
    % Parameters:
    %  - pos  :     distribution of sites
    %  - sam  :     initial sample
    %  - nsim :     number of changes before convergence
    %  - rngb :     radius of neighbourhood
    %  - rocc :     radius of region occupied by a site
    % Given an initial condition, at each MC step we want t converge to a
    % more likely state. To do so, from x, we move 1 receptor (k) at a 
    % time to one of its neigbouring states(k+1) with probability 
    % P = min(Nk+1/Nk,1) where Nk is the amount of neighbours that the kth
    % site has.
    
    N = length(pos(1,:));
    M = length(sam(1,:));
    
    ngb_n = zeros(1,N);
    sam_i = find(ismember(pos',sam','rows')); % indices of the positions that are sampled
    
    
    for jj = 1:N
            d = dist(pos(:,jj)',pos);
            ingb = find(d<rngb & d>0);
            ngb_n(jj) = length(ingb);
            % save neighbours for each site, and their corresponding indeces
            ngb_v{jj} = [pos(:,ingb);ingb];
    end
    
    
    for j=1:nsim
        for k = 1:M
            psam = sam(:,k); % point sampled
            kk   = sam_i(k);
            ngb_h = ngb_v{kk};
            if ~isempty(ngb_h)
                chc      = randsample(ngb_n(kk),1);
                ngb_cdt = ngb_h(3,chc); % site candidate for moving the receptor (returns the index)
                % check for overlap
                if M~=1 
                    d   = dist(pos(:,ngb_cdt)',pos(:,sam_i(sam_i~=kk)));
                    m   = min(d);
                else m = rocc+1;
                end
                if m>rocc %(no overlap)
                    p = min(1,ngb_n(kk)/ngb_n(ngb_cdt));
                    rdm = rand();
                    if rdm < p
                        sam_i(k) = ngb_cdt;
                    end
                end
            end
        end
    end
    
    y = pos(:,sam_i);
    
end