% Script for creating plots of growth coupling in supplementary Fig 3

solver='gurobi';
modelName='e_coli_core.mat';
saveData=0;
targetRxn='EX_succ_e';


% initialise COBRA toolbox
changeCobraSolver(solver,'all',0);
model=readCbModel(modelName);

% find indices of target/biomass
targetInd=find(ismember(model.rxns,targetRxn));
growthInd=find(model.c==1);

% find gene-associated reactions
geneAssocList=find(~ismember(model.rules,''));

% create LP model and store as a parallel pool constant
lpModel=convertModel(model);
lpModel.c(targetInd)=-10^-4;
poolModel=parallel.pool.Constant(lpModel);

% test all single deletions
startDels=false(length(geneAssocList));
for a=1:length(geneAssocList)
    startDels(a,a)=true;
end
startFitness = myFitness3obj(startDels, poolModel, geneAssocList, targetInd, growthInd, length(model.rxns), 10^-8, 0, 10^-5, -inf, 10^-4, false, [], []);
startFitness = - startFitness;

% find essential reactions and remove them from further consideration- no
% point in testing combinations of essential reactions
essRxns=startFitness(:,3)==-inf;
nonEssList=geneAssocList(~essRxns);
nonEssFitness=startFitness(~essRxns,:);

% test all combinations of non-essential reactions
prodVals=zeros(length(nonEssList));
csVals=zeros(length(nonEssList));



for a=1:length(nonEssList)
    disp(a/length(nonEssList))
    
    % store metrics determined in startFitness
    prodVals(a,a)=nonEssFitness(a,2);
    csVals(a,a)=nonEssFitness(a,3);
    
    % test combinations of deletion a with all other target deletions
    if a==length(nonEssList)
        % skip if a is the final deletion- no combinations left to test
        continue
    end
    
    % create deletion matrix
    testDels=false(length(nonEssList)-a,length(nonEssList));
    testDels(:,a)=true;
    for b=1:size(testDels,1)
        testDels(b,a+b)=true;
    end
    
    % test deletions
    testFitness = myFitness3obj(testDels, poolModel, nonEssList, targetInd, growthInd, length(model.rxns), 10^-8, 0, 10^-5, -inf, 10^-4, false, [], []);
    testFitness = - testFitness;
    
    % store result in matrix
    prodVals(a,a+1:end)=testFitness(:,2);
    csVals(a,a+1:end)=testFitness(:,3);
    
    % Also store result in opposite direction (Deletion A + Deletion B == Deletion B + Deletion A)
    prodVals(a+1:end,a)=testFitness(:,2)';
    csVals(a+1:end,a)=testFitness(:,3)';
    
end

% sort deletions alphabetically

rxnList=model.rxns(nonEssList);
[rxnList,sortInd]=sort(rxnList);

prodVals=prodVals(sortInd,sortInd);
csVals=csVals(sortInd,sortInd);

% show deletions that cause coupling
[xCoupled,yCoupled]=find(prodVals>0);
t=table(rxnList(xCoupled),rxnList(yCoupled),zeros(length(xCoupled),1),zeros(length(xCoupled),1),'VariableNames',{'Rxn1','Rxn2','Prod','CouplingStrength'});
for a=1:size(t,1)
    t{a,3}=prodVals(xCoupled(a),yCoupled(a));
    t{a,4}=csVals(xCoupled(a),yCoupled(a));
end
t=sortrows(t,3);
disp(t)

% set custom colormap

colPoints=[0,0.01,1;
    0.55,0.8,0;
    0.3,0.9,0.6];

colPoints2=[colPoints(1,1):(colPoints(1,2)-colPoints(1,1))/30:colPoints(1,2), colPoints(1,2):(colPoints(1,3)-colPoints(1,2))/30:colPoints(1,3);
    colPoints(2,1):(colPoints(2,2)-colPoints(2,1))/30:colPoints(2,2), colPoints(2,2):(colPoints(2,3)-colPoints(2,2))/30:colPoints(2,3);
    colPoints(3,1):(colPoints(3,2)-colPoints(3,1))/30:colPoints(3,2), colPoints(3,2):(colPoints(3,3)-colPoints(3,2))/30:colPoints(3,3);]';

colPoints3=[linspace(0,1,size(colPoints2,1))',linspace(0,1,size(colPoints2,1))',linspace(0,1,size(colPoints2,1))'];

% Plot product synthesis

colormap([colPoints3;colPoints2])
% adjust zeros so they appear white instead of coloured
scaledProd=prodVals;
scaledProd(scaledProd==0)=-0.0001;

figure()
heatmap(rxnList,rxnList,scaledProd,'colormap',[colPoints3;colPoints2],'ColorLimits',[-max(max(scaledProd)),max(max(scaledProd))]);
set(gcf,'color','w')
set(gca,'title','Product synthesis')
set(gca,'FontSize',7)


% Plot coupling strength

% Start normalising coupling strength to max/min values so
% positive/negative values are correct colours
scaledCs=csVals;

% Dividing positive coupling strength by max possible coupling strength
scaledCs(scaledCs>0)=scaledCs(scaledCs>0)/2; 

% Dividing negative coupling strengths by lowest value identified
tempVals=unique(csVals(:));
if tempVals(1)~=-inf
    minCoupling=tempVals(1);
else
    minCoupling=tempVals(2);
    scaledCs(scaledCs<minCoupling)=minCoupling;
end
scaledCs(scaledCs<0)=scaledCs(scaledCs<0)/-minCoupling;

% Making values with zero coupling strength appear white rather than coloured (as not growth
% coupled)
scaledCs(scaledCs==0)=-0.0001;

% Plot normalised values
figure()
heatmap(rxnList,rxnList,scaledCs,'colormap',[colPoints3;colPoints2],'ColorLimits',[-1,1]);
set(gcf,'color','w')
set(gca,'title','Coupling Strength')
set(gca,'FontSize',7)

