function x = sample_min_dist(pos, n, dis)
    % Arguments
    %   - pos = set of possible positions (2-by-N)
    %   - n   = sample size
    %   - dis = minimal distance between sample points
    % This function chooses a sample of [n] elements from the possible
    % positions [pos], with the constraint that chosen points are at least
    % a distance [dis] away from each other.
    % The function returns a 2-by-n matrix with the positions of sampled
    % points as column vectors, or an empty matrix, if it is impossible to
    % find an appropriate distribution of sampled points.
    
    N = length(pos(1,:));
    
    in = 1:N;    % vector of possible indeces
    sam_i = randsample(N,n);         % choose n of the N possibilities
    sam = pos(:,sam_i);
    com_i = in(~ismember(in,sam_i));      % indeces of non-chosen indeces
    com = pos(:,com_i);
    
    % Check if sample points are far enough from one another
    i=2;
    fail = 0;
    tr_max = 500;
    while i<=n && fail == 0
        n_tr = 1;
        d = dist(sam(:,i)',sam(:,1:i-1));
        min_d = min(d);
        
        while min_d < dis && n_tr < tr_max
            temp = sam_i(i);
            new = randsample(N-n,1);
            sam_i(i) = com_i(new);
            com_i(new) = temp;
            sam(:,i) = pos(:,sam_i(i));
            com(:,new) = pos(:,temp);
            
            d = dist(sam(:,i)',sam(:,1:i-1));
            min_d = min(d);
            n_tr = n_tr+1;
        end
        if n_tr == tr_max
            %disp('Could not find appropriate distribution.')
            fail = 1;
            x=[];
        end 
        i=i+1;
    end
    
    if fail == 0
        x = sam;
    end    
end