function ex = examineGcOptResults(results, results_wt)
% Examine gcOpt results struct for GC mutants
% Maximal growth rate
% Minimal production rate at maximal growth (Ymin)
% Ratio between pRmin and maximal theoretical production rate
% Minimal productivity at maximal growth
% Ratio between Pmin and maximal theoretical productivity
% GCS
% Growth state according to MiMBl

% if gcOpt results are invalid stop examination
if results.validity ~= 1
    ex.gcs   = -2;
    return
end

rxnNumSubs              = find(ismember(results.model.rxns,results.model.subsRxn));
rxnNumTarget            = find(ismember(results.model.rxns,results.model.targetRxn));

yRange_wt      = results_wt.yieldR;
pRange_wt      = results_wt.prodR;

yRange_mut       = results.printData.yieldR;
mRange_mut       = results.printData.muR;
pRange_mut       = results.printData.prodR;



%% Wild type
[th_maxY, ~]    = max(yRange_wt);
[th_maxP, ~]    = max(pRange_wt);

compSave    = [];
%% Maximal growth rate
[ex.maxMu, ~]   = max(mRange_mut);
compSave           = [compSave,ex.maxMu];

%% Minimal production rate at maximal growth (minPR)
isMax               = yRange_mut(mRange_mut==ex.maxMu);
ex.minPR_maxMu    = min(isMax);
compSave            = [compSave,ex.minPR_maxMu];

%% Ratio between PRmin and maximal theoretical production rate [%]
ex.ratio_minPR   = (ex.minPR_maxMu/th_maxY)*100;
compSave           = [compSave,ex.ratio_minPR];

%% Minimal productivity at maximal growth (Pmin)
isMax           = pRange_mut(mRange_mut==ex.maxMu);
ex.minP_maxMu   = min(isMax);
compSave           = [compSave,ex.minP_maxMu];

%% Ratio between Pmin and maximal theoretical productivity [%]
ex.ratio_minP   = (ex.minP_maxMu/th_maxP)*100;
compSave           = [compSave,ex.ratio_minP];

%% Growth coupling strength
ex.gcs  = calcGCS(results, results_wt);
compSave           = [compSave,ex.gcs];

% compact structure
ex.compSave     = compSave;

%% Yield at maximal growth
model   = results.model;
% apply KOs
model   = changeRxnBounds(model,model.rxns(results.KORxnNum(:,1)),0,'b');
% FBA
model   = changeObjective(model,model.bmRxn);
sol     = optimizeCbModel(model,'max','one');
ex.Y_maxMu   = sol.x(rxnNumTarget)/(-sol.x(rxnNumSubs));



end