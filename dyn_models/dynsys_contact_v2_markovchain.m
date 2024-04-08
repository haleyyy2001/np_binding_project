function y = dynsys_contact_v2_markovchain(kon,koff,v,nt,pia)
    % Given the parameters of the NP/TCR cluster system, this function
    % returns the equilibrium distribution of the states as well as the
    % mean residency times for each sate.
    % The parameters of the system are as follows
    %  - kf : forward rate in a single, mononeric TCR/pMHC interaction
    %  - kb : backward rate in a single, monomeric TCR/pMHC interaction
    %  - v : valence of the NP
    %   vh is the effective (/contact) valence of the considered NP
    %  - nt : number of TCR in the specific contact area
    vh = round(v/7);
    N = min(vh,nt);
    
    % 
    if N>1
        if N>length(pia)
            pia = [pia,zeros(1,N-length(pia))];
        end
        
        Pia = pia(2:N); % (the first entry for pia is always 1)
    
    
        % Compute the forward rate for each state
        Pia;
        f = vh*Pia.*(nt*ones(1,N-1)-(1:N-1))*kon;
        f0 = v*kon;
        f = [f,0];
        % Compute the backward rate for each state
        b = (1:N)*koff;
    
        % Compute forward and bckward rates as transition probabilities
        F = f./(f+b);
        B = b./(f+b);
        
        matss = zeros((N+1)/2*(2*nt-N+2));
        %mats = cell(1,nt+1);
        
        for k =1
            mats = [0,1];
            cord = (k-1)*k/2+1;
            [li,co] = size(mats);
            matss(cord:cord+li-1,cord:cord+co-1) = mats;
        end
        
        for k = 2:N
            D1 = ((0:k-2)./(nt-(k-1:-1:1)).*F(k-1:-1:1))';
            D1 = [D1;0];
            D2 = (B(k-1:-1:1))';
            D2 = [D2;0];
            D3 = ((nt-(k-1))./(nt-(k-1:-1:1)).*F(k-1:-1:1))';
            D3 = [D3;0];
            De = [zeros(k-1,1);1];
            mats = spdiags([D1,D2,D3,De],[1,-1,-k,0],2*k-1,k)';
            cord = (k-1)*k/2+1;
            [li,co] = size(mats);
            matss(cord:cord+li-1,cord:cord+co-1) = mats;
        end
        
        for k = N+1:nt
            D1 = ((k-N-1:k-2)./(nt-(N:-1:1)).*F(N:-1:1))';
            D1 = [D1;0];
            D2 = (B(N:-1:1))';
            D2 = [D2;0];
            D3 = ((nt-(k-1))./(nt-(N:-1:1)).*F(N:-1:1))';
            D3 = [D3;0];
            De = [zeros(N,1);1];
            mats = spdiags([D1,D2,D3,De],[1,-1,-N,0],2*N,N+1)';
            cord = (N+1)*N/2+1+(k-N-1)*(N+1);
            [li,co] = size(mats);
            matss(cord:cord+li-1,cord:cord+co-1) = mats;
        end
        
        for k = nt+1
            D1 = F(N:-1:1)';
            D1 = [D1;0];
            D2 = (B(N:-1:1))';
            D2 = [D2;0];
            De = [zeros(N,1);1];
            mats = spdiags([D1,D2,De],[1,-1,0],N+1,N+1)';
            cord = length(matss(:,1))-N;
            [li,co] = size(mats);
            matss(cord:cord+li-1,cord:cord+co-1) = mats;
        end
        
        matQ = matss((diag(matss)~=1),(diag(matss)~=1));
        matemp = sparse(eye(length(matQ(:,1))) - matQ);
        matN = inv(matemp);
%         for i = 1:length(matQ(:,1))
%             tp = [zeros(i-1,1);1;zeros(length(matQ(:,1))-i,1)];
%             matN(:,i) = matemp\tp;
%         end
        
        matR = matss((diag(matss)~=1),(diag(matss)==1));
        
        matB = matN(1,:)*matR;
        
        matssbi = [[matss((diag(matss)~=1),(diag(matss)~=1)),...
            matss((diag(matss)~=1),(diag(matss)==1))];...
            [matss((diag(matss)==1),(diag(matss)~=1)),...
            matss((diag(matss)==1),(diag(matss)==1))]];
        
%         lhs_ss = matss - eye((N+1)/2*(2*nt-N))
%         lhs_ss = sparse(lhs_ss);
%         rhs_ss = zeros((N+1)/2*(2*nt-N),1);
%         ss = lhs_ss\rhs_ss;
%         % Include both pieces of information in the same structure
%         y = ss;

        rtdc = b+f;
        rtdu = -[0,f(1:end-1)];
        rtdl = -[b(2:end),0];
    
        matrt = spdiags([rtdl',rtdc',rtdu'],[-1,0,1],N,N);
        rhs_rt = ones(N,1);
    
        rt = matrt\rhs_rt;
    
        if abs(sum(matB(1,:))-1)<=0.001
            ss = matB(1,:)';
        else
            ss = zeros(size(matB(1,:)));
            ss(end) = 1;
            ss = ss(:);
        end
        
        
        y = struct('ss',ss,'rt',rt(1));
        
    else
        f0 = v*kon;
        % Compute the backward rate for each state
        b = koff;
    
        ss = 1;
        rt = 1/b;
    
        y = struct('ss',ss,'rt',rt(1));
    end
end
    
    