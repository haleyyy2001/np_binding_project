function [y, ma, md] = louis_activationEngage(kd,KD,theta_deg,np_params,ifng_params,pho_params,pp)
%% ------------------------------------------------------------------------
% Initializing Parameters
% -------------------------------------------------------------------------
    M = 1;%
    Mtot=5000;
	koff=25;            % Was 10 before
    KOFF = 25;
    kon = koff / kd;
    KON = KOFF / KD;
% Experimental parameters .................................................
    np_params = num2cell(np_params);
    [r,v,d] = np_params{:};
    
    theta = theta_deg * pi / 180;
    rc = r * sin(theta);                %Radius of contact area on Tcell surface
    SAeff = 2*pi*(r^2)*(1-cos(theta));  %Effective surface area of NP 
    SAtot = 4*pi*(r^2);                 %Total surface area of NP
    Acont = pi * r ^ 2;                 %Contact area on Tcell Surface
    vh = round(SAeff * v / SAtot);      %Effective valency of NP
    
% Phosphorylation Parameters
   % pho_params = num2cell(pho_params);
   % [phi,b,St,C_star,Npho] = pho_params{:};
   pt = pho_params;
    
% IFNg Activation Parameters
    ifng_params = num2cell(ifng_params);
    [mu,alph,C] = ifng_params{:};
    
%% ------------------------------------------------------------------------
% Load Data
%--------------------------------------------------------------------------

    %pp = ['insProb_results_Manuela/'];
    %pp = ['louis_temp/theta',num2str(theta_deg),'/'];
    %pp = ['louis_results/dTCR10/'];
pp = ['louis_results/'];
% Load simulation data for NP .............................................
    load([pp,'insprob_r',num2str(r),'_v',num2str(v),'.mat'])
    load([pp,'insprobNP_r',num2str(r),'.mat'])
    
% Determine the carrying capacities of nano-clusters ......................
    Mvec = M*probBindNPi(:);
    
%% ------------------------------------------------------------------------
%  Contact Area Steady State
% -------------------------------------------------------------------------

% 1. Find the maximum number of TCR that can be bound by a single NP.......
    maxTCRcov = 1;
    for i = 1:length(distCov)
        maxTCRcov = max(maxTCRcov,length(distCov{i}));
    end
    
% 2. Compute the unbinding time of NPs and the average number of bound TCR per NP
    avgs_contAreaModel = dynsys_means_v2_markovchain(koff/kd,koff,v,maxTCRcov,insProbTCRtoNP);
    backwardRate = unbindingRate(distCov,avgs_contAreaModel(:,3));
    
%% ------------------------------------------------------------------------
% Find the steady state from the nano-cluster model
% -------------------------------------------------------------------------
    stStateNC = louis_stState_ncModel(koff/kd,koff,KON,KOFF,backwardRate,Mvec,r,v,d,[pp,'Continuation_Results/']);
    %stStateNC = louis_stState_ncModel(koff/kd,koff,backwardRate,Mvec,r,v,d,'zeroNewton/'); % For Manuela's results
                                                                                           %Replace pp with 'zeroNewton/' for Manuela's Continuation Results.
    stStateNC =  stStateNC/sum(stStateNC);
    
%% ------------------------------------------------------------------------
% Cooperativity in the nano-cluster model
   %{
    f = zeros(1, length(Mvec));
    for i=1:length(Mvec)
        f(i) = v * KON * kon * d * (Mvec(i)-stStateNC(i+1));
        b(i) = KOFF * backwardRate(i);
    end
    
    disp(stStateNC)
    
    disp(log(f ./ b) / log(f(1)/b(1)))
    disp(length(Mvec))
    
    alphaN = cumsum(log(f ./ b)) ./ [1:length(Mvec)];  % Cooperativity of NC model with N = Mvec nanoparticles.
    disp(alphaN)
    alphaN = alphaN(end) / log(f(1)/b(1));
    
    disp(alphaN)
 
    
    figure(1)
    plot([1:length(Mvec)],f)
    ylabel('Forward rate (f)')
    xlabel('Number of NPs in NC')
    
    figure(2)
    plot([1:length(Mvec)],b)
    xlabel('Number of NPs in NC')
    ylabel('Backward rates (b)')
    
    figure(3)
    plot([1:length(Mvec)], cumsum(log(f ./ b)) ./ [1:length(Mvec)])
    title('Average association constants of n-valent interaction')
    ylabel('Log(K)')
    hline(0)
    xlabel('Number of NPs in NC')
    
    figure(4)
    plot([1:length(Mvec)], cumsum(log(f./b)) ./ [1:length(Mvec)] / log(f(1)/b(1)));
    hline(0)
    title('Cooperativity of NC model')
    ylabel('Alpha')
    xlabel('Number of NPs in NC')
    
    figure(5)
    plot([0:length(Mvec)], stStateNC)
    title('Occupancy of NC model steady states')
    xlabel('Number of NPs in NC')
%}
%% ------------------------------------------------------------------------
% Finding average number of bound TCR
% -------------------------------------------------------------------------

    % Determine the average number of bound TCR per TCRnc based on the
    % order of binding (i.e. #bound by first, second etc)
    bdNumber = [];
    for i = 1:length(distCov)
        bdNumber = [bdNumber, distCov{i}*avgs_contAreaModel(1:length(distCov{i}),2)];
    end
    % Sum to find total bound number per NC with j NPs attached
    bdNumber = cumsum(bdNumber);
    
    % Average bound per NC
    convAvg = stStateNC(2:end)'*bdNumber(:)/sum(stStateNC(2:end));
    
%% ------------------------------------------------------------------------
% Fitting Distributions (Optional)
% -------------------------------------------------------------------------
    %{
    if (round(vh)>1 & length(avgs_contAreaModel(:,2))>1)
        disp(length(avgs_contAreaModel(:,2)))
        [alphas, beta]  = dist_fit_gamma(distCov,avgs_contAreaModel(:,2),'Gamma');
    
        alphas = cumsum(alphas);
        alpha = alphas'*Mtot*stStateNC(2:end);
        
        meanBoundTCR = alpha * beta;         % Mean bound TCR
        sdBoundTCR = sqrt(alpha * beta^2);   % Standard Deviation of bound TCR
        
        maxConvValues = meanBoundTCR + 5 * sdBoundTCR;            % Max and Min used to generate linspace
        minConvValues = max(meanBoundTCR - 5 * sdBoundTCR, 0);
        stpConvValues = (maxConvValues - minConvValues)/200;
        convValues = [linspace(0,minConvValues,200),linspace(minConvValues,maxConvValues,200)];
        convDensit = pdf('gamma',convValues,alpha,beta);           % Generate pdf distribution of bound TCRs
        
        ma = convAvg;                           % Mean bound TCR per NC from NC steady state model
        md = meanBoundTCR/Mtot;                 % Mean bound TCR per NC from distribution fitting
        
        
        trigTCR = convValues*pt;                            % Converting from bound TCRs to triggered TCRs
        convActiv = ifngProd(trigTCR,mu,alph,C);            % Converting from triggered TCRs to IFNg production
        activlev = convDensit*convActiv'*stpConvValues;     % Integrating IFNg production over pdf distribution
        
        %{
        activlev = convDensit * convValues' * stpConvValues;                % Number of Bound TCRs
        activlev = louis_Phosphorylation(activlev,phi,b,St,C_star,koff,Npho);  % Number of Triggered TCRs
        activlev = ifngProd(activlev,mu,alph,C);                            % IFNg Production
        %}
    else
        activlev = convAvg*Mtot*(1-stStateNC(1));                           % Number of bound TCRs
        %activlev = louis_Phosphorylation(activlev,phi,b,St,C_star,koff,Npho);  % Number of triggered TCRs
        activlev = pt*activlev;
        activlev = ifngProd(activlev,mu,alph,C);                            % Take (1-stStateNC(1) for prob of bound NP
        ma = convAvg;
        md = convAvg;
        
    end
    %}
    activlev = convAvg*Mtot*(1-stStateNC(1));
    activlev = pt*activlev;%lev number of bound tcr
    activlev = ifngProd(activlev,mu,alph,C);
    
    ma = convAvg;
    md = convAvg;
    
    y = activlev/ 1000;
end
    
    
    
    