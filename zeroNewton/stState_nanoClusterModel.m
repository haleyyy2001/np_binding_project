function y = stState_nanoClusterModel(kon,koff,B,Mvec,r,v,d)
    % ,ntcr,dtcr)
    % kd is the dissociation constant
    % B is the backward rate
    % Mvec is the NP capacity of the nano-cluster
    % r is the radius of the nanoparticle
    % v is the valence of the NP
    % d id the dose of NP
    % Need input arguments for ntcr and dtcr
    
    kd = koff/kon;
    R = kon*v*d;
    
%     load(['contKdSteady_r',num2str(r),...
%         'v',num2str(v),'d',num2str(round(1000*d)),'ntcr',...
%         num2str(ntcr),'dtcr',num2str(round(1000*dtcr)),'.mat'])
    load(['contKdSteady_r',num2str(r),...
        'v',num2str(v),'d',num2str(round(1000*d)),'.mat'])
    
    [clKd,inKd] = min(abs(kdGrid/kd-1));
    statePrevious = stateKdContGrid(:,inKd);
    if length(statePrevious)>length(Mvec)
        Mvec = [Mvec;zeros(length(statePrevious)-length(Mvec)-1,1)];
        B = [B;ones(length(statePrevious)-length(B)-1,1)];
    end
    dX = 1      ;
    while norm(dX) > 10^-6
        dX = newtonZeros(R,Mvec,B,statePrevious);
        stateNow = statePrevious+dX;
        statePrevious = stateNow;
    end
   y = stateNow; 
end