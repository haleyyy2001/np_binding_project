function [fitDists_singleA, fitDists_singleB]  = dist_fit_gamma(dists,stStatesTCR,distop)
%% Setup
posNP = length(dists);
%distop = 'Gamma';
%distop = 'Normal';
stStatesTCR = stStatesTCR(:);

%% Computing the cumulative distributions
startPoints = zeros(1,2);
for i = 1:posNP
    cumDists{i} = cumsum(dists{i});
    
    if i==1
    lD = length(dists{i});
    mu = stStatesTCR(1:lD)'*dists{i}(:);                     %Causing issues at Angle=30, values were 1 and 0 anyway.
    vv = (stStatesTCR(1:lD).^2)'*dists{i}(:) - mu^2;
    
    startPoints(i,1) = mu^2/vv;
    startPoints(i,2) = vv/mu;
    end
end
%% Determining the abscisse for the cumulative distributions
absc1 = [(stStatesTCR(2:end)+stStatesTCR(1:end-1))/2;...
    stStatesTCR(end)+(stStatesTCR(end)-stStatesTCR(end-1))/2];

%% Fit the cumulative distribution function with a Hill function
ftfct = @(a,b,x)cdf(distop,x,a,b);
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,Inf],...
               'StartPoint',[1 1]);
ft = fittype(@(a,b,x)ftfct(a,b,x),'independent', 'x','options',fo);
fitAll1 = [];
startPoint = startPoints;
for i = 1:posNP
    fo = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0,0],...
                   'Upper',[Inf,Inf],...
                   'StartPoint',startPoint);%startPoints(:,i)');
    ft = fittype(@(a,b,x)ftfct(a,b,x),'independent', 'x','options',fo);
    lD = length(cumDists{i});
    temp = fit([0;absc1(1:lD)],[0,cumDists{i}]',ft);
    fitAll1 = [fitAll1,[temp.a; temp.b]];
    startPoint = [temp.a; temp.b];
    
end
    
    %['m(beta) = ', num2str(mean(fitAll1(2,:))),'; t1 = ',num2str(t1)]
%%
if distop == 'Gamma';
    allB = fitAll1(2,:);
    mnB = mean(allB);
    varB = var(allB);
    %['Variance in beta is ',num2str(sqrt(varB)/mnB*100),'%.']
    
    hiB = mnB+1.96*sqrt(varB);
    loB = mnB-1.96*sqrt(varB);
    
    fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[Inf],...
               'StartPoint',[1]);
    ftfct = @(a,x)cdf(distop,x,a,mnB);
    ft = fittype(@(a,x)ftfct(a,x),'independent', 'x','options',fo);
    fitAllB = [];
    
    %ftfct = @(a,x)cdf(distop,x,a,mnB);
    %fitAllB = [];
    for i = 1:posNP
        
        fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[Inf],...
               'StartPoint',[fitAll1(1,i)]);
        ft = fittype(@(a,x)ftfct(a,x),'independent', 'x','options',fo);
        
        lD = length(cumDists{i});
        temp = fit([0;absc1(1:lD)],[0,cumDists{i}]',ft);
        fitAllB = [fitAllB;temp.a];
    end
    
end
%% Results
%toc
fitDists_singleA = fitAllB;
fitDists_singleB = mnB;
end

