function results_m = init_gcOpt(model,probOpts,genOpts)
% Initialize solving MILP to find product-growth-coupled strain designs
%
% INPUT
%     model:              Stoichiometric metabolic model in COBRA format
%                             Additional necessary fields:
%                             - bmRxn:        Identifier of biomass formation reaction
%                             - subsRxn:      Identifier of substrate uptake reaction (exchange reaction)
%                             - targetRxn     Identifier of product formation reaction
%                             
%     probOpts:           Options for gcOpt
%                             Necessary fields:
%                             - maxKO:    maximal number of deletions
%                             Optional fields:
%                             - fixMu:        objective growth rate
%                             - notKORxns:    Identifiers of reactions excluded from target space
%
%     genOpts:            General (technical) options for gcOpt
%                             - compressFlag:         (1) Compress model by removing blocked reactions
%                                                     (0) No compression (default)
%                             - TimeLimit:            Solver time limit in seconds (default: Inf)
%                             - printFlag:            (0) Do not plot results (default)
%                                                     (1) Plot results as a production envelope
%                             - solvThreads:          Number of threads to be utilized
%                             - solvCuts:             gurobi solver cut setting (default: 1) 
%                             - Presolve              gurobi solver presolve settings (default: 1)
%                             - avBnd:                auxiliary variable bounds (default: 500)                    
%  
% OUTPUT
%     results_m:          results structure
%                             - KORxnNum:           reaction numbers of deletion targets in the original model
%                             - KORxnNames:       reaction names of deletion targets
%                             - KOs:              deletion targets assigend in vector format
%                             - fluxes:           flux distribution at optimality of the MILP
%                             - rawRes:           unprocessed results
%                             - wildtypeData:     production envelope data of wildtype model
%                             - printData:        production envelope data of GC mutant
%                             - model:            metabolic model in COBRA format
%                             - ex:               key data of mutant model
%                                                 .maxMu:         Maximal growth rate
%                                                 .minPR_maxMu:   minimally guaranteed target production rate at maximal growth rate
%                                                 .ratio_minPR:   ratio between minPR_maxMu and theoretical maximal target production rate
%                                                 .minP_maxMu:    minimally guaranteed productivity (mmol/h) at maximal growth
%                                                 .ratio_minP:    ratio between minP_maxMu and maximal theoretical productivity
%                                                 .gcs:           grwoth-coupling strength
%                                                 .Y_maxMu:       target product yield at maximal growth rate
% 
% AUTHOR:
%     - Tobias Alter      29th January 2018
%         Institute of Applied Microbiology, RWTH Aachen University

%%

fprintf(['\n\t\t\t  GGG    CCC',...
         '\n\t\t\t G      C      OOO    PPP   TTTTT',...
         '\n\t\t\t GGGG   C     O   O  P   P    T',...
         '\n\t\t\t G   G  C     O   O  P   P    T',...
         '\n\t\t\t  GGG    CCC   OOO   PPPP     T',...
         '\n\t\t\t                     P',...
         '\n\t\t\t                     P\n'])


%% Check specified options
if nargin < 2
    probOpts    = [];
    genOpts     = [];
end
if nargin < 3
    genOpts     = [];
end
if isempty(probOpts)
    probOpts    = [];
end

%% Check problem related options
if ~isfield(probOpts,'optType')
    % 0: gcOpt
    probOpts.optType    = 0;
    optType             = probOpts.optType;
else
%     probOpts.optType    = 0;
    optType     = probOpts.optType;
end

if ~isfield(probOpts,'sense')
    probOpts.sense    = 1;
else
    if strcmp(probOpts.sense,'max')
        probOpts.sense  = 1;
    elseif strcmp(probOpts.sense,'min')
        probOpts.sense  = -1;
    else
        error('Unknown assignment of optimization direction (probOpts.sense)')
    end
end


if ~isfield(probOpts,'maxKO')
    error('No intervention set size specified. Set probOpts.maxKO!')
else
    maxKO       = probOpts.maxKO;
end

if ~isfield(probOpts,'fixMu')
    probOpts.fixMu       = [];
end


if ~isfield(probOpts,'bmRxn')
    error('Specify biomass-reaction-identifier!')
else
    bmRxn           = probOpts.bmRxn;
    model.bmRxn     = bmRxn;
end

if ~isfield(probOpts,'subsRxn')
    error('Specify substrate-uptake-reaction-identifier (C-source)!')
else
    subsRxn         = probOpts.subsRxn;
    model.subsRxn   = subsRxn;
end

if ~isfield(probOpts,'targetRxn')
    error('Specify objective-reaction-identifier!')
else
    targetRxn       = probOpts.targetRxn;
    model.targetRxn     = targetRxn;
end

if ~isfield(probOpts,'notKORxns')
    notKORxns   = [];
else
    notKORxns       = probOpts.notKORxns;
end


%% Check general options
if ~isfield(genOpts,'compressFlag')
    reductionFlag   = 0; 
    % 0: No reduction (not available)
    % 1: Remove blocked reactions (not available)
else
    reductionFlag   = genOpts.compressFlag;
end

if ~isfield(genOpts,'solver')
    genOpts.solver.LP       = 0;
    genOpts.solver.MILP     = 0;
    genOpts.solver.MINLP    = 0;
    % 0: Gurobi Solver (MILP)
else
    if ~isfield(genOpts.solver,'LP')
        genOpts.solver.LP   = 0;
    end
    if ~isfield(genOpts.solver,'MILP')
        genOpts.solver.MILP   = 0;
    end
    if ~isfield(genOpts.solver,'MINLP')
        genOpts.solver.MINLP   = 0;
    end
end
solver  = genOpts.solver;

if ~isfield(genOpts,'modelType')
    genOpts.modelType       = 0;
else
    if isempty(genOpts.modelType)
        genOpts.modelType       = 0;
    end
end

if ~isfield(genOpts,'printFlag')
    genOpts.printFlag       = 1;
end

if ~isfield(genOpts,'solvThreads')
    genOpts.solvThreads       = 0;
end

if ~isfield(genOpts,'solvCuts')
    genOpts.solvCuts       = 1;
end

if ~isfield(genOpts,'avBnd')
    genOpts.avBnd       = 500;
end

if ~isfield(genOpts,'TimeLimit')
    genOpts.TimeLimit       = Inf;
end

if ~isfield(genOpts,'Presolve')
    genOpts.Presolve       = 1;
end


%% Preprocess Model
% Reduce maximal upper and lower boundaries to 1000 or -1000, respectively
model.ub(model.ub>1000)     = 1000;
model.lb(model.lb<-1000)    = -1000;

%% Handle parallelization
p   = gcp('nocreate');
if isempty(p)
    if genOpts.solvThreads <= 0
        % start default parallel cluster
        parpool;
    else
        % start paralel cluster
        parpool(genOpts.solvThreads);
    end
else
    if (p.NumWorkers ~= genOpts.solvThreads) && (genOpts.solvThreads>0)
        delete(p)   % close parallel cluster
        parpool(genOpts.solvThreads)    % start new parallel cluster
    end
end
p                       = gcp('nocreate');
genOpts.solvThreads     = p.NumWorkers;
disp(['Number of parallel workers: ',num2str(p.NumWorkers)])



%% Optimisation algorithm initialisation
% Model to MIP
[A_j, b_j, lb_j, ub_j,cO, vSize, zSize, lSize, ySize, B_flux, mapRxnsRed, B, nonKnock, model_s,results_wt] =  ...
    model2MIP(model,notKORxns,bmRxn,targetRxn,subsRxn,maxKO,reductionFlag,solver,probOpts,genOpts);


% Solve Problem
 [results] = solveGcOptMILP(model_s, A_j,b_j,lb_j,ub_j,vSize, zSize, lSize, ySize,B,B_flux,...
                            mapRxnsRed,nonKnock,results_wt,solver,genOpts);

                        
% Postprocessing
results_m = matchCompr2OrigMod(results,model_s,model,probOpts);  

results_m.rawRes          = results;
results_m.rawRes.model_s  = model_s;
results_m.wildtypeData    = results_wt;
% Test KO strategy
results_m =   testStrategy(model, results_m);



%% Print solution
% Yield over Mu
if ~isempty(results_m.KOs)
    printData              = YoM(model,bmRxn,subsRxn,targetRxn,results_m.KORxnNum,genOpts.printFlag);
end

results_m.printData     = printData;
results_m.model         = model;

%% Examine gcOpt Result
fprintf('--> Examine results \n')
if (optType==0) && isfield(results_m.wildtypeData,'yieldR') && isfield(results_m.printData,'yieldR') 
    results_m.ex = examineGcOptResults(results_m, results_m.wildtypeData);
end



