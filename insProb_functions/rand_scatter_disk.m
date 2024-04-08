function x = rand_scatter_disk(r, n, dis)
    % Arguments:
    %   - r:        radius of the disk on which to sample,
    %   - n:        # of points to be sampled
    %   - dis:      minimum distance required between sampled points
    % This function generates a distribuion of [n] points on a disk of
    % radius [r]. The sampled points need to be at least a distance [dis]
    % away from each other.
    % This function returns the distribution of sampled points in a 2-by-n array.
    % If it is impossible to find a distribution that satisfies the
    % constraints on the integral, the function returns an empty array.
    
    tr_max = 500;
    
    th = 2*pi*rand(1,n);                        %Random angle
    rho = rand(1,n)+rand(1,n);                  %Random radius
    rho(rho>1) = 2-rho(rho>1);
    x = [r*rho.*cos(th);r*rho.*sin(th)];
    
    
    if nargin > 2
    i=2;    
    fail=0;
        while i <= n && fail == 0

            n_tr = 1;
            d     = dist(x(:,i)',x(:,1:i-1));
            min_d   = min(d);
            while min_d < dis && n_tr < tr_max
                rho = rand(1,1)+rand(1,1);
                rho(rho>1) = 2-rho(rho>1);
                th = 2*pi*rand(1,1);
                x(:,i) = [rho*r*cos(th),rho*r*sin(th)];
                d     = dist(x(:,i)',x(:,1:i-1));
                min_d   = min(d);
                n_tr = n_tr+1;
            end
            if n_tr == tr_max
                %disp('Could not find available position.')
                fail=1;
                x = [];
            end
            i=i+1;
        end
        
    end
end
