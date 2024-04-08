function y = activationFunction(r,v,kd,kp,kappa)

    %% ------------------------------------------------------------------------
    % Initialization of parameters
    % -------------------------------------------------------------------------
    
    % Physiological parameters ................................................
%     kd = .01;
%     koff = 1;
%     kp = 10;
%     kappa=7;
%     M = 1;%
    Mtot=20000;
    koff=1;
%     % Experimental parameters .................................................
%     r = 14;
%     v = 31;
%     d = 1.875;
    %% ------------------------------------------------------------------------
    % Load stored data
    %--------------------------------------------------------------------------
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Make this into a function that can
    % be used at each step. Cleans code
    % and likely speeds the code
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    % Load simulation data for NP .............................................
    load(['insprob_r',num2str(r),'_v',num2str(v),'.mat'])
    load(['insprobNP_r',num2str(r),'.mat'])
    
    % Determine the carrying capacities of nano-clusters ......................
    Mvec = M*probBindNPi(:);
    
    %% ------------------------------------------------------------------------
    % Contact area model steady state
    % -------------------------------------------------------------------------
    
    % 1. Find the maximum number of TCR that can be bound by a single NP.......
    maxTCRcov = 1;
    for i = 1:length(distCov)
        maxTCRcov = max(maxTCRcov,length(distCov{i}));
    end
    
    % 2. Compute the unbinding time and the average number of bound TCR per NP
    %       orig syntax: dynsys_means(kon,koff,v,maxCov,insProbTCRtoNP)........
    %       Since now know only dependent on koff and Kd, can write:...........
    avgs_contAreaModel = dynsys_means(koff/kd,koff,v,maxTCRcov,insProbTCRtoNP);
    backwardRate = unbindingRate(distCov,avgs_contAreaModel(:,3));
    
    %% ------------------------------------------------------------------------
    % Find the steady state from the nano-cluster model
    % -------------------------------------------------------------------------
    stStateNC = stState_nanoClusterModel(koff/kd,koff,backwardRate,Mvec,r,v,d);
    
    %% ------------------------------------------------------------------------
    % Distribution fitting (distribution of the number of bound TCR per T-cell)
    % -------------------------------------------------------------------------
    [alphas, beta]  = dist_fit_gamma(distCov,avgs_contAreaModel(:,2));
    alpha = alphas'*Mtot*stStateNC(2:end);
    
    mnBoundTCR = alpha*beta;
    sdBoundTCR = sqrt(alpha*beta^2);
    maxConvValues = mnBoundTCR+3*sdBoundTCR;
    minConvValues = max(mnBoundTCR-3*sdBoundTCR,0);
    stpConvValues = (maxConvValues - minConvValues)/200;
    convValues = linspace(minConvValues,maxConvValues,stpConvValues);
    convDensit = pdf('gamma',convValues,alpha,beta);
    %% ------------------------------------------------------------------------
    % Activation
    % -------------------------------------------------------------------------
    trigTCR = convValues*(kp/(kp+koff))^kappa;
    convActiv = ifngProd(trigTCR);
    activlev = convDensit*convActiv'*stpConvValues;
    
    
    y = activlev;
end







