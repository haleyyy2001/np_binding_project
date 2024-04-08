function y = fitPiMuAlpha(chain,sschain)
    
    
    srchain = sort(sschain);
    maxer = srchain(round(length(sschain)/5));
    
    newchain = chain(sschain<= maxer,:);
    
    fto = fitoptions('Method','nonlinearleastsquares','lower',[eps eps eps],'upper',[inf inf inf],'startpoint',[20 3 1])
    ftp = fittype(@(a,b,c,x1,x2) a*((x1.^b).*(x2.^c)).^(-1),'independent',{'x1','x2'},'options',fto)
    [fre, gof] = fit([newchain(:,2),newchain(:,3)],newchain(:,4),ftp);
    
    plt = plot(fre,[newchain(:,2),newchain(:,3)],newchain(:,4));
    set(plt(1),'edgecolor','none','facealpha',0.5);
    set(plt(2),'markerfacecolor',[255 153 153]/255);
   
    xlabel('\pi')
    ylabel('m')
    zlabel('\alpha')
    rmse = 0;
    
    title({'\alpha = a \times (\pi^b m^c)^{-1}',...
        ['a = ',num2str(fre.a),', b = ',num2str(fre.b),', c = ',num2str(fre.c)],...
        ['RMSE = ',num2str(gof.rmse), ', r = ', num2str(gof.rsquare)]})
    
    y = fre;
end