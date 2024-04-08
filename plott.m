
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
