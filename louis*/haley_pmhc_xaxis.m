
data=get_ifndata(14)
data_54=data(1:7,2:3);
data_31=data(8:11,2:3);
data_14=data(12:17,2:3);
data_11=data(18:23,2:3);
data_8=data(24:28,2:3);

%plot(data_54(:,1),data_54(:,2))
%hold on
%plot(data_31(:,1),data_31(:,2))
%plot(data_14(:,1),data_14(:,2))
%plot(data_11(:,1),data_11(:,2))
%plot(data_8(:,1),data_8(:,2))
%hold off
%Hill(x,y,noParam,maximum,slope,halfActiv,intercept);
%Hill(data_54(:,1), data_54(:,2),3,3,10,2.5,0)
x_54=54*(data_54(:,1))
x_31=31*(data_31(:,1))
x_14=14*(data_14(:,1))
x_11=11*(data_11(:,1))
x_8=8*(data_8(:,1))
% MAPPING: Emax = b(1),  EC50 = b(2)
err_ee=[];
for qc = 1:10
    hill_fit = @(b,x)  b(1).*(x.^qc)./((b(2))^qc+(x.^qc));
    b0 = [4; 9];                                  % Initial Parameter Estimates
    B = lsqcurvefit(hill_fit, b0, x_54 ,data_54(:,2));
    AgVct = linspace(min(x_54), max(x_54 ));   % Plot Finer Resolution
    B_31 = lsqcurvefit(hill_fit, b0, x_31 ,data_31(:,2));
    AgVct_31 = linspace(min(x_31 ), max(x_31 ));   % Plot Finer Resolution
    B_14 = lsqcurvefit(hill_fit, b0, x_14 ,data_14(:,2));
    AgVct_14 = linspace(min(x_14 ), max(x_14 ));
    B_11 = lsqcurvefit(hill_fit, b0, x_11 ,data_11(:,2));
    AgVct_11 = linspace(min(x_11 ), max(x_11 ));
    B_8 = lsqcurvefit(hill_fit, b0,  x_8 ,data_8(:,2));
    AgVct_8 = linspace(min(x_8) , max(x_8));
    di_54=data_54(:,2)-hill_fit(B,x_54);
     di_31=data_31(:,2)-hill_fit(B_31,x_31);
      di_14=data_14(:,2)-hill_fit(B_14,x_14);
       di_11=data_11(:,2)-hill_fit(B_11,x_11);
        di_8=data_8(:,2)-hill_fit(B_8,x_8);

    err_ee(qc)=norm(di_54)+norm(di_31)+norm(di_14)+norm(di_11)+norm(di_8);

end
hill_fit = @(b,x)  b(1).*(x.^2)./((b(2))^2+(x.^2));
    b0 = [4; 9];                                  % Initial Parameter Estimates
    B = lsqcurvefit(hill_fit, b0, x_54 ,data_54(:,2));
    AgVct = linspace(min(x_54), max(x_54 ));   % Plot Finer Resolution
    B_31 = lsqcurvefit(hill_fit, b0, x_31 ,data_31(:,2));
    AgVct_31 = linspace(min(x_31 ), max(x_31 ));   % Plot Finer Resolution
    B_14 = lsqcurvefit(hill_fit, b0, x_14 ,data_14(:,2));
    AgVct_14 = linspace(min(x_14 ), max(x_14 ));
    B_11 = lsqcurvefit(hill_fit, b0, x_11 ,data_11(:,2));
    AgVct_11 = linspace(min(x_11 ), max(x_11 ));
    B_8 = lsqcurvefit(hill_fit, b0,  x_8 ,data_8(:,2));
    AgVct_8 = linspace(min(x_8) , max(x_8));
figure()
    plot(AgVct, hill_fit(B,AgVct),'-r')
hold on

legend('54')
plot(x_54 ,data_54(:,2), 'rp','HandleVisibility','off')
plot(x_31 ,data_31(:,2), 'bp','HandleVisibility','off')
plot(x_14(:,1),data_14(:,2), 'yp','HandleVisibility','off')
plot(x_11(:,1),data_11(:,2), 'gp','HandleVisibility','off')
plot(x_8(:,1),data_8(:,2), 'mp','HandleVisibility','off')
plot(AgVct_31, hill_fit(B_31,AgVct_31), '-b','DisplayName','31')

plot(AgVct_14, hill_fit(B_14,AgVct_14), '-y','DisplayName','14')

plot(AgVct_11, hill_fit(B_11,AgVct_11), '-g','DisplayName','11')

plot(AgVct_8, hill_fit(B_8,AgVct_8), '-m','DisplayName','8' )
 
grid
xlabel('pmhc')
ylabel('infg')

hold off
err_ee
text=["ec50","54:",B(2),"31:",B_31(2),"14:",B_14(2),"11:",B_11(2),"8:",B_8(2)]


