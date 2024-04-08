function y = activationFunction_v2_markovchain(r,v,d,kd,pi,mu,alph,C)

    %% ------------------------------------------------------------------------
    % Initialization of parameters
    % -------------------------------------------------------------------------
    
    % Physiological parameters ................................................
%     kd = .01;
%     koff = 1;
%     kp = 10;
%     kappa=7;
    M = 1;%
   Mtot=5000;
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
    avgs_contAreaModel = dynsys_means_v2_markovchain(koff/kd,koff,v,maxTCRcov,insProbTCRtoNP);
    backwardRate = unbindingRate(distCov,avgs_contAreaModel(:,3));
    
    %% ------------------------------------------------------------------------
    % Find the steady state from the nano-cluster model
    % -------------------------------------------------------------------------
    stStateNC = stState_nanoClusterModel(koff/kd,koff,backwardRate,Mvec,r,v,d);
    
    %% ------------------------------------------------------------------------
    % Determine the average number of bound
    % -------------------------------------------------------------------------
    bdNumber = [];
    
    for i = 1:length(distCov)
        bdNumber = [bdNumber, distCov{i}*avgs_contAreaModel(1:length(distCov{i}),2)];
    end
    stStateNC;
    bdNumber = cumsum(bdNumber);
    
    convValues = stStateNC(2:end)'*bdNumber(:)/sum(stStateNC(2:end));
    
    %% ------------------------------------------------------------------------
    % Distribution fitting (distribution of the number of bound TCR per T-cell)
    % -------------------------------------------------------------------------
%     [alphas, beta]  = dist_fit_gamma(distCov,avgs_contAreaModel(:,2));
%     alpha = alphas'*Mtot*stStateNC(2:end);
%     
%     mnBoundTCR = alpha*beta;
%     sdBoundTCR = sqrt(alpha*beta^2);
%     maxConvValues = mnBoundTCR+3*sdBoundTCR;
%     minConvValues = max(mnBoundTCR-3*sdBoundTCR,0);
%     stpConvValues = (maxConvValues - minConvValues)/200;
%     convValues = linspace(minConvValues,maxConvValues,200);
%     convDensit = pdf('gamma',convValues,alpha,beta);
    
%     %% ------------------------------------------------------------------------
%     % Figure to check steady states
%     % -------------------------------------------------------------------------
%     vecPlot = (1:20000);
%     if nargin == 5
%         actiPlot = ifngProd(vecPlot);
%     elseif nargin == 6
%         actiPlot = ifngProd(vecPlot,mu);
%     elseif nargin == 8 
%         actiPlot = ifngProd(vecPlot,mu,alph,C);
%     end
%     figure(25)
%     plot(vecPlot,actiPlot/max(actiPlot),'linewidth',2);
%     hold on
%     distPlot = gampdf(vecPlot,alpha,beta);
%     plot(vecPlot*pi,distPlot/max(distPlot),'linewidth',2)
%     hold off
%     pause(1)
    %% ------------------------------------------------------------------------
    % Activation
    % -------------------------------------------------------------------------
    
    
    trigTCR = convValues*pi*5000;
    if nargin == 5
        convActiv = ifngProd(trigTCR);
    elseif nargin == 6
        convActiv = ifngProd(trigTCR,mu);
    elseif nargin == 8 
        convActiv = ifngProd(trigTCR,mu,alph,C);
    end
    activlev = 25000*convActiv;
    
    
    y = activlev/1000;%*25000/1000;
end







