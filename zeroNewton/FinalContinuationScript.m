clear all
close all

r14 = get_ifndata(14);
r20 = get_ifndata(20);

rpos = [14*ones(length(r14(:,1)),1);20*ones(length(r20(:,1)),1)];
vpos = [r14(:,1);r20(:,1)];
dpos = [r14(:,2);r20(:,2)];
k=1;
for i = 1:length(rpos)%8:11;%
    r = rpos(i);
    v = vpos(i);
    D = dpos(i);
    %% Initialize parameters involved in continuation
    load(['insprob_r',num2str(r),'_v',num2str(v),'.mat'])      % Insertion probability of TCR onto NP
    load(['insprobNP_r',num2str(r),'.mat'])                    % Insertion probability of NP onto nano-cluster
    M = 1;%20000;                                              % Total number of nano-clusters
    distcum = probBindNPi(:);                                  % Distribution of the probability of reaching a given state
    Mvec = M*distcum;  %                                       % Prop of nano-clusters that can reach a given state
    kk=1;
    %% Continuation parameters
    konInitial = .0001;                   % Initial value for on-rate (known to integrate well)
    koffInitial = 0.1;               % Initial value for off-rate (known to integrate well)
    kdInitial = koffInitial/konInitial;
    kdMin = 1e-2;
    kdMax = 1e4;
    contStep = 1.15;
    newtonTol = 1e-6;
    %% Find steady state for initial (kon,koff) pair
    temp = dynsys_means(konInitial,koffInitial, v, 10,insProbTCRtoNP);
    B = unbindingRate(distCov,temp(:,3));
    
    % Solve an ODE to find the initial point from which to start the
    % continuation
    R = konInitial*v*D;
    opts = odeset('Events',@(x,t)eventfun(R,Mvec,B,x,t));
    [t,y] = ode15s(@(t,x)f(R,Mvec,B,x,t),[0,500],[M,zeros(1,length(insProbNPtoNC))]');%,opts);
    initState = y(end,:);
    
    
    % Run continuation with respect to parameter kd
    [kdGrid, stateKdContGrid] = continuation(3,kdMin,kdMax,contStep,konInitial,koffInitial,...
        initState,newtonTol,D,v,insProbTCRtoNP,distCov,Mvec);
    
    
    cols = parula(12);
    
    %if v ==31
        
    figure(1);
    %subplot(2,2,k)
    
    for ii = 1 :length(initState)
        semilogx(kdGrid,20000*stateKdContGrid(ii,:),'linewidth',1,'color',cols(ii,:))
        hold on
    end
    grid on
    axis([-inf inf 0 20000*M])
    xlabel('K_D')
    ylabel('st. state')
    title(['d = ',num2str(D)])
    pause(.01)
    hold off    
    k=k+1;
    if k==5
        legend('Y_0','Y_1','Y_2','Y_3','Y_4','Y_5','Y_6','Y_7','Y_8','Y_9')
        saveas(figure(1),['contKD_example'],'epsc')
    end
    %end
    pause(0.5)
    save(['contKdSteady_r',num2str(r),'v',num2str(v),'d',num2str(1000*D)],'kdGrid','stateKdContGrid')
end