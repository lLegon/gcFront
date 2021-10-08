%% function to test viability og gcOpt KO strategy
function results = testStrategy(model, results)
targetRxnNum    = find(ismember(model.rxns,model.targetRxn));

% set up FBA model
model   = changeRxnBounds(model,model.bmRxn,results.wildtypeData.optBmFlux,'b');
model   = changeRxnBounds(model,model.rxns(results.KORxnNum(:,1)),0,'b');
% FBA
model   = changeObjective(model, model.targetRxn);
sol     = optimizeCbModel(model,'min');
if strcmp(sol.origStat,'OPTIMAL')
    targetRxn   = sol.x(targetRxnNum);
    if abs(targetRxn) < 0.0001
        % KO strategy insufficient
        val = 0;
    else
        val = 1;
    end
else
    val = 0;
end
results.validity   = val;   


end