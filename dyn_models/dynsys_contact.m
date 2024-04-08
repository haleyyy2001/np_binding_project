function y = dynsys_contact(kon,koff,v,vh,nt,pia)
    % Given the parameters of the NP/TCR cluster system, this function
    % returns the equilibrium distribution of the states as well as the
    % mean residency times for each sate.
    % The parameters of the system are as follows
    %  - kf : forward rate in a single, monomeric TCR/pMHC interaction
    %  - kb : backward rate in a single, monomeric TCR/pMHC interaction
    %  - v : valence of the NP
    %   vh is the effective (/contact) valence of the considered NP
    %  - nt : number of TCR in the specific contact area
    N = min(vh,nt);
    
    if N>1
        if N>length(pia)
            pia = [pia,zeros(1,N-length(pia))];
        end
        
        Pia = pia(2:N); % (the first entry for pia is always 1)
    
    
        % Compute the forward rate for each state
        Pia;
        f = vh*Pia.*(nt*ones(1,N-1)-(1:N-1))*kon;
        f0 = v*kon*nt;                                                      % Originally use v, but testing vh for first CAM transition.
        f = [f0,f];
        % Compute the backward rate for each state
        b = (1:N)*koff;
    
        % Computation of steady distribution
        %   Generate matrix to find steady state
        ssdc = [-f,0]+[0,-b];                   % center diagonal for SS
        ssdu = [0,b];                           % upper diagonal for SS
        ssdl = [f,0];                           % lower diagonal for SS
    
        %   Here, don't use spdiags since bottom line is all ones
        matss = spdiags([ssdu',ssdc',ssdl'],[1,0,-1],N+1,N+1);
        matss(end,:) = ones(1,N+1);
    
        %   Solve for steady state distribution. The first (N-1) equations
        %   represent the computation of dXi/dt = 0, while the last one
        %   corresponds to the equation X0 + X1 + ... + XN = 1
    
        rhs_ss = [zeros(N,1);1];
        ss = matss\rhs_ss;
    
        % Computation of mean residency times
        rtdc = b+[f(2:end),0];
        rtdu = -[0,f(2:end)];
        rtdl = -[b(2:end),0];
    
        matrt = spdiags([rtdl',rtdc',rtdu'],[-1,0,1],N,N);
        rhs_rt = ones(N,1);
    
        rt = matrt\rhs_rt;
    
    
        % Include both pieces of information in the same structure
        y = struct('ss',ss,'rt',rt);
    else
        f = v*kon;                                                          % Originally use v, replaced with vh to test
        b = koff;
        ssdc = [-f,0]+[0,-b];                   % center diagonal for SS
        ssdu = [0,b];                           % upper diagonal for SS
        ssdl = [f,0];                           % lower diagonal for SS
        matss = spdiags([ssdu',ssdc',ssdl'],[1,0,-1],N+1,N+1);
        matss(end,:) = ones(1,N+1);
        rhs_ss = [zeros(N,1);1];
        ss = matss\rhs_ss;
        
        rtdc = b+[f(2:end),0];
        rtdu = -[0,f(2:end)];
        rtdl = -[b(2:end),0];
    
        matrt = spdiags([rtdl',rtdc',rtdu'],[-1,0,1],N,N);
        rhs_rt = ones(N,1);
    
        rt = matrt\rhs_rt;
        
        y = struct('ss',ss,'rt',rt);
    
    end
end
    
    