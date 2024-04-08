%% Pre-processing (lol)

data = [[repmat(14,28,1),get_ifndata(14)];[repmat(20,21,1),get_ifndata(20)]];

all_sims = cell(length(data(:,1)),4);

for i = 1:length(data(:,1))
    r = data(i,1);
    v = data(i,2);
    d = data(i,3);
    load(['insprobNP_r',num2str(r),'.mat']);
    load(['insprob_r',num2str(r),'_v',num2str(v),'.mat']);
    load(['contKdSteady_r',num2str(r),...
        'v',num2str(v),'d',num2str(1000*d),'.mat']);
    
    all_sims{i,1} = insProbTCRtoNP;
    all_sims{i,2} = distCov;
    all_sims{i,3} = probBindNPi;
    all_sims{i,4} = stateKdContGrid;
end