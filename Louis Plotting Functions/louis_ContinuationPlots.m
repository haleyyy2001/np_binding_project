% Continuation Methods Plots to determine steady states of the TNM and
% their deendence on the parameter Kd.

pp = ['louis_results/Continuation_Results/'];
r = 14;
v=14;
d=1.5;                                                      % Same dose for both v-types.
NTCR = [12,17,22,27];

f1 = figure(1);
i=1;

for ntcr = NTCR
    load([pp,'contKdSteady_r',num2str(r),'v',num2str(v),'d',num2str(d*1000),'.mat'])%'ntcr',num2str(ntcr),

    cols = parula([length(stateKdContGrid(:,1))]);                              %Colors for plot

    a=subplot(2,2,i);

    for ii = 1 :length(stateKdContGrid(:,1))
        txt = ['$\hat{Y}_{', num2str(ii-1),'}$'];
        semilogx(kdGrid,20000*stateKdContGrid(ii,:),'linewidth',2,'color',cols(ii,:),'DisplayName',txt)
        hold on
    end
    hold off
    grid on
    
    axis([-inf inf -1000 20000])
    a.YAxis.Exponent=4;
    xlabel('K_D')
    ylabel('Steady State')
    title(['v = ', num2str(v),', nTCR = ',num2str(ntcr),])
    set(gca,'FontSize',16,'LineWidth',2)
    
    i=i+1;
end

leg = legend('Position',[0.45 0.1 1 1],'Interpreter','latex','Fontsize',20);
title(leg, 'Bound NPs')

ax1 = subplot(2,2,1);
ax2 = subplot(2,2,3);

AddLetters2Plots({ax1,ax2},{'A','B'},'Hshift',-0.05,'Vshift',-0.05,'FontSize',32) 

f1.Units='centimeters';
f1.Position=[1 50 50 30];
