function results = YoM(model,bmRxn,subsRxn,objRxn,KORxnNum,plotFlag)

% close all


results = [];

% Objective reaction number
objRxnNum       = find(strcmp(model.rxns,objRxn));
% objRxn          = model.rxns(objRxnNum);
bmRxnNum        = find(strcmp(model.rxns,bmRxn));
subsRxnNum      = find(strcmp(model.rxns,subsRxn));

flags    = 0;
defModelLb                  = model.lb;
defModelUb                  = model.ub;


for i=1:length(KORxnNum)
    
    defModelLb                  = model.lb;
    defModelUb                  = model.ub;
    % Set Knockouts
    model.lb(KORxnNum(i))  = 0;       
    model.ub(KORxnNum(i))  = 0;

    % % Set Knockdowns (bounds)
    % if any(results.KDs)
    %     KDPos       = find(results.KDs);
    %     for i=1:size(KDPos,1)
    %         pos             = KDPos(i);
    %         model.lb(pos)     = model.lb(pos)*results.KDs(pos);
    %         model.ub(pos)     = model.ub(pos)*results.KDs(pos);
    %     end
    % end

    %% Check model feasibility
    % Set objective (max production rate with ensured growth)
    model = changeRxnBounds(model,bmRxn,0,'l'); %set growth rate minimum
    model           = changeObjective(model,objRxn);
    % Conduct FBA
    FBAsol          = optimizeCbModel(model,'max');

    if strcmp(FBAsol.origStat,'INFEASIBLE')
        flags    = 1;
        disp('Infeasible intervention strategy')
        model.lb    = defModelLb;
        model.ub    = defModelUb;
    else
        defModelLb                  = model.lb;
        defModelUb                  = model.ub;
        flags    = 0;
%         break;
    end
end

if flags
    error('No feasible intervention strategies!')
end


%% Define yield range
% fixSubs     = -0.5;
% model = changeRxnBounds(model,subsRxn,fixSubs,'b'); %fix substrate uptake rate
model   = changeRxnBounds(model,bmRxn,0,'l'); %set growth rate minimum
model   = changeObjective(model,model.rxns(objRxnNum));
FBAsol  = optimizeCbModel(model,'max');

maxVIntervent           = FBAsol.x(objRxnNum);


if isempty(maxVIntervent)
    warning('No production with selected intervention strategy!')
    return
end

yieldVecInterventLoop       = [0:0.02:1]'.*maxVIntervent;
muRangeIntervent            = zeros(size(yieldVecInterventLoop));
prodIntervent               = muRangeIntervent;
yieldVecIntervent           = muRangeIntervent;


%% Calculate yield over growth rate

model           = changeObjective(model,bmRxn);

% Alternative plotting routing
% % Loop through yield values
% for i=1:size(yieldVecInterventLoop,1)
% 
%     % lock yield
%     model = changeRxnBounds(model,objRxn,yieldVecInterventLoop(i),'b'); %anaerobic
% 
%     % Conduct FBA
%     FBAsol          = optimizeCbModel(model,'max');
% 
%     if ~strcmp(FBAsol.origStat,'INFEASIBLE')
%         % Calculate yield
%         yieldVecIntervent(i)         = FBAsol.x(objRxnNum)/(-1*FBAsol.x(subsRxnNum));
%         muRangeIntervent(i)          = FBAsol.x(bmRxnNum);
%         subsUptakeIntervent(i)       = FBAsol.x(subsRxnNum);
%         % calculate productivity
%         prodIntervent(i)             = FBAsol.x(objRxnNum)*(exp(muRangeIntervent(i))-1);
%     else
%         yieldVecIntervent(i)    = -1;
%         muRangeIntervent(i)     = -1;
%         prodIntervent(i)        = -1;
%     end
% 
% end


%% Different approach maximise/minimise growth rate

model.lb    = defModelLb;
model.ub    = defModelUb;

fixSubs     = model.lb(subsRxnNum);
% % fixSubs     = -0.5;
% model           = changeObjective(model,bmRxn);
% FBAsol          = optimizeCbModel(model,'max');
% fixSubs         = FBAsol.x(subsRxnNum)

model = changeRxnBounds(model,subsRxn,fixSubs,'b'); %fix substrate uptake rate

model           = changeObjective(model,bmRxn);
% Conduct FBA
FBAsol          = optimizeCbModel(model,'max');

if strcmp(FBAsol.origStat,'INFEASIBLE')
    model = changeRxnBounds(model,subsRxn,0,'u');
    FBAsol          = optimizeCbModel(model,'max');
    if strcmp(FBAsol.origStat,'INFEASIBLE')
        error('Infeasible model')
    end
end

muRmin      = [0:0.02:1]'.*FBAsol.x(bmRxnNum);
muRmax      = muRmin;
yieldRmin   = zeros(size(muRmin));
yieldRmax   = yieldRmin;
prodRmin    = yieldRmin;
prodRmax    = yieldRmin;
fluxMin     = yieldRmin;
fluxMax     = yieldRmin;
bmYieldMin     = yieldRmin;
bmYieldMax     = yieldRmin;

model           = changeObjective(model,objRxn);
for i=1:size(muRmax,1)

    % lock yield
    model = changeRxnBounds(model,bmRxn,muRmin(i),'b'); 

    % Conduct FBA
    FBAsolMax          = optimizeCbModel(model,'max');
    FBAsolMin          = optimizeCbModel(model,'min'); 

    if ~strcmp(FBAsolMax.origStat,'INFEASIBLE')
        fluxMax(i)              = FBAsolMax.x(objRxnNum);
        % Calculate yield
        yieldRmax(i)            = FBAsolMax.x(objRxnNum)/(-1*FBAsolMax.x(subsRxnNum));
        bmYieldMax(i)           = muRmax(i)/(-1*FBAsolMax.x(subsRxnNum));
        subsUptakeIntervent(i)       = FBAsolMax.x(subsRxnNum);
        % calculate productivity
        prodRmax(i)             = FBAsolMax.x(objRxnNum)*(exp(muRmax(i))-1);
    else
        yieldRmax(i)    = -1;
        prodRmax(i)     = -1;
        fluxMax(i)      = -1;
        bmYieldMax(i)   = -1;
    end
    
    if ~strcmp(FBAsolMin.origStat,'INFEASIBLE')
        fluxMin(i)              = FBAsolMin.x(objRxnNum);
        % Calculate yield
        yieldRmin(i)            = FBAsolMin.x(objRxnNum)/(-1*FBAsolMin.x(subsRxnNum));
        bmYieldMin(i)           = muRmin(i)/(-1*FBAsolMin.x(subsRxnNum));
        subsUptakeIntervent(i)  = FBAsolMin.x(subsRxnNum);
        % calculate productivity
        prodRmin(i)             = FBAsolMin.x(objRxnNum)*(exp(muRmin(i))-1);
    else
        yieldRmin(i)    = -1;
        prodRmin(i)     = -1;
        fluxMin(i)      = -1;
        bmYieldMin(i)   = -1;
    end

end

swapVec     = linspace(size(muRmax,1),1,size(muRmax,1));

yieldR      = [yieldRmin;yieldRmax(swapVec)];
prodR       = [prodRmin;prodRmax(swapVec)];
muR         = [muRmin;muRmax(swapVec)];
flux        = [fluxMin;fluxMax(swapVec)];
bmYield     = [bmYieldMin;bmYieldMax(swapVec)];

% transform units
bmYield     = bmYield.*1000;    % g/mol

% Delete infeasible solutions
del     = (yieldR<0)|(prodR<0)|(muR<0);

yieldR(del)     = [];
prodR(del)      = [];
muR(del)        = [];
flux(del)       = [];
bmYield(del)    = [];

% [yieldR,I]  = sort(yieldR,'ascend');
% prodR       = prodR(I);
% muR         = muR(I);


if plotFlag
    fprintf('--> Print results \n')
   %     printYieldoverMu(muR,yieldR,prodR,model.rxnNames{objRxnNum},0);
    figure;
%     printYieldoverMu_add(muR,yieldR,prodR,model.rxnNames{objRxnNum});
     printYieldoverMu_add(bmYield,yieldR,prodR,model.rxnNames{objRxnNum});

end


% figure
% hold on
% printYieldoverMu_add(muRmin,yieldRmin,prodRmin, [objRxn,' max/min mu']);
% printYieldoverMu_add(muRmax,yieldRmax,prodRmax, [objRxn,' max/min mu']);
% hold off

%% Plot results
% % Delete infeasible solutions
% yieldVecIntervent(muRangeIntervent<0)    = [];
% prodIntervent(muRangeIntervent<0)        = [];
% muRangeIntervent(muRangeIntervent<0)     = [];
% 
% [minYield,idx]      = min(yieldVecIntervent);
% if minYield~=0
%     disp('Minimum Yield is not zero')
%     addMu               = [0:(muRangeIntervent(idx)/10):muRangeIntervent(idx)]';
%     yieldVecIntervent   = [yieldVecIntervent;ones(size(addMu,1),1).*minYield];
%     muRangeIntervent    = [muRangeIntervent;addMu];
%     prodIntervent       = [prodIntervent;(-subsUptakeIntervent(idx)*minYield).*(exp(addMu)-1)];
% end

% if plotFlag
%     printYieldoverMu(muRangeIntervent,yieldVecIntervent,prodIntervent, objRxn,0);
% end


results.yieldR  = yieldR;
results.prodR   = prodR;
results.muR     = muR;
results.flux    = flux;
results.bmYield = bmYield;
results.Mu_Y_P_F_BMY  = [muR,yieldR,prodR,flux,bmYield];








end