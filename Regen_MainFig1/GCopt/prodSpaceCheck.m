function output = prodSpaceCheck(model,bmRxn,subsRxn,objRxn,sense)
% Explore unperturbed (productivity) solution space

fprintf('\t ...Determine Solution Space Edges \n')

bmRxnNum        = find(strcmp(model.rxns,bmRxn));
subsRxnNum      = find(strcmp(model.rxns,subsRxn));
objRxnNum       = find(strcmp(model.rxns,objRxn));

% Optimization direction
if sense == 1
    FBAsense    = 'max';
elseif sense == -1
    FBAsense    = 'min';
end

%% Check model feasibility
% Set objective 
model = changeRxnBounds(model,bmRxn,0.01,'l'); %set growth rate minimum
model           = changeObjective(model,bmRxn);
% Conduct FBA
FBAsol          = optimizeCbModel(model,'max');

if strcmp(FBAsol.origStat,'INFEASIBLE')
    save error
    error('Infeasible model')
end


% %% Define yield range
% % fixSubs     = -10;
% % model = changeRxnBounds(model,subsRxn,fixSubs,'b'); %fix substrate uptake rate
% model = changeRxnBounds(model,bmRxn,0,'l'); %set growth rate minimum
% model           = changeObjective(model,objRxn);
% FBAsol          = optimizeCbModel(model,FBAsense);
% 
% maxVIntervent           = FBAsol.x(objRxnNum);
% 
% yieldVecLoop       = [0:0.02:1]'.*maxVIntervent;
% muRange            = zeros(size(yieldVecLoop));
% prod               = muRange;
% yieldVec           = muRange;
% subsUptake         = muRange;
% objFlux            = muRange;
% %% Calculate yield over growth rate
% 
% model           = changeObjective(model,bmRxn);
% 
% % Loop through yield values
% for i=1:size(yieldVecLoop,1)
% 
%     % lock yield
%     model = changeRxnBounds(model,objRxn,yieldVecLoop(i),'b'); %anaerobic
% 
%     % Conduct FBA
%     FBAsol          = optimizeCbModel(model,FBAsense);
% 
%     if ~strcmp(FBAsol.origStat,'INFEASIBLE')
%         objRxnFlux          = abs(FBAsol.x(objRxnNum));
%         % Calculate yield
%         yieldVec(i)         = objRxnFlux;%/(-1*FBAsol.x(subsRxnNum));
%         muRange(i)          = FBAsol.x(bmRxnNum);
%         subsUptake(i)       = FBAsol.x(subsRxnNum);
%         % calculate productivity
%         prod(i)             = objRxnFlux*(exp(muRange(i))-1);
%         % save objective flux rate
%         objFlux(i)          = objRxnFlux;
%     else
%         yieldVec(i)    = -1;
%         muRange(i)     = -1;
%         prod(i)        = -1;
%         subsUptake(i)  = -1;
%     end
% 
% end
% 
% % Delete infeasible solutions
% yieldVec(yieldVec==-1)    = [];
% muRange(muRange==-1)      = [];
% prod(prod==-1)            = [];
% subsUptake(subsUptake==-1)= [];





% output.yieldR           = yieldVec;
% output.muR              = muRange;
% output.prodR            = prod;
% output.subsUpR          = subsUptake;
% output.objFluxR         = objFlux;


results                 = YoM(model,bmRxn,subsRxn,objRxn,[],0);
output.yieldR           = results.yieldR;
output.muR              = results.muR;
output.prodR            = results.prodR;
output.bmYield          = results.bmYield;
output.Mu_Y_P_F_BMY     = results.Mu_Y_P_F_BMY;
% output.subsUpR          = subsUptake;
% output.objFluxR         = objFlux;


end