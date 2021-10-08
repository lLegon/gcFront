function [] = testFastPros
%testFastPros is a test program contains script that test FastPros function.
% 
% Aug. 5th, 2013    Satoshi OHNO

%% Input

% Metabolic model to be used
load('Ecoli_core_model.mat','model')

% Target reaction to be maximized
targetRxn = 'EX_succ(e)';

% Biomass reaction as objective function
biomassRxn = 'Biomass_Ecoli_core_w_GAM';

% Exchange raction of oxygen uptake
oxygenRxn = 'EX_o2(e)';


%% Optional input

% Load results of FVA or not in model reduction
options.loadFVAFlux = false;

% Reaction list as knockout candidates
switch model.description
    case {'Ecoli_core_model'}
        options.rxnList = setdiff(model.rxns,...
            {'ENO','GAPD','PGK','PGM',...       % Essential reaction in actual cell.
            'H2Ot','EX_h2o(e)',...
            'EX_acald(e)','ACALDt',...
            'EX_co2(e)','CO2t'});               % Hard to knockout only these reaction.
end

% Maximum knockout number in FastPros
options.maxKoNum =5;

% Knockout strain number to be selected as parent strain in each iteration
options.selStrainNum = [];

% Select only strains whose utarget increased by the current knokcout
options.selIncUtargetStrains = true;

% Verbose flag
options.verbFlag = false;
% options.verbFlag = true;

if options.verbFlag == true
    disp( 'Model:        E.coli core model' )
    disp(['TargetRxn:    '  targetRxn])
    disp(['BiomasssRxn:  '  biomassRxn])
    disp(['OxygenRxn:    '  oxygenRxn])
    disp(['MaxKoNum:     '  num2str(options.maxKoNum)])
    disp(' ')
end

%% Other inputs required for this test progaram

% Load reduced model or not
loadReducedModel = false;
% loadReducedModel = true;

% Load FastPros result or not
loadFPResult = false;



%% Change reaction bounds of the model

% Maximum uptake rate of oxygen molecule
OUR = 5;
% Exchange reaction of carbon soruce uptake
carbonSourceInfo.rxn = 'EX_glc(e)';
% Maximum uptake rate of carbon soruce
carbonSourceInfo.uptakeRate = 10;
% ATP maintenance reaction
atpMaintenaceRxn = 'ATPM';
% ATP maintenance cost
NGAM = 8.39;
% Minimun growth rate as threshhold
minGrowth = 0.05;

model = changeRxnBounds(model,carbonSourceInfo.rxn,-carbonSourceInfo.uptakeRate,'l');
model = changeRxnBounds(model,oxygenRxn,-OUR,'l');
model = changeRxnBounds(model,atpMaintenaceRxn,NGAM,'l');
model = changeRxnBounds(model,biomassRxn,minGrowth,'l');
model = changeObjective(model,{biomassRxn,targetRxn},[1,10^-5]);
if strcmp(model.description,'Ecoli_core_model')
    %Change objective function to maximize biomass production
    model.c = zeros(length(model.rxns),1);
    biomassRxnID = findRxnIDs(model,biomassRxn);
    model.c(biomassRxnID) = 1;
end

%% Reduce the model for FastPros

modelOri = model;
if loadReducedModel == 1;
    load(['reduced_' model.description])
    biomassRxn = model.rxns{1};
    targetRxn = model.rxns{2};
    oxygenRxn = model.rxns{3};
else
    [model,biomassRxn,targetRxn,oxygenRxn] = ...
        reduceModelForFP(modelOri,biomassRxn,targetRxn,oxygenRxn,options);
    save(['reduced_' model.description],'model')
end
if options.verbFlag == true
    disp(' ')
end


%% FastPros

if loadFPResult == 1
    load('testFastProsResult')
else
    [FastProsSolution, FastProsStatus, FastProsResult] = ...
        FastPros(model, biomassRxn, targetRxn, oxygenRxn,options);
%     save('testFastProsResult', 'FastProsSolution', 'FastProsStatus', 'FastProsResult')
end


%% Check FastPros result

loadData = load('testFastProsResult');

isFailure = 1;
if ~isempty('FastProsSolution')
    if length(FastProsSolution) == length(loadData.FastProsSolution)
        isFailure = 0;
        for i = 1 : length(FastProsSolution)
            if length(FastProsSolution(i).prodRates) ~= length(loadData.FastProsSolution(i).prodRates)
                isFailure = isFailure+1;
                break
            end
            if any(FastProsSolution(i).prodRates ~= loadData.FastProsSolution(i).prodRates)
                isFailure = isFailure+1;
            end
            if any(FastProsSolution(i).utarget ~= loadData.FastProsSolution(i).utarget)
                isFailure = isFailure+1;
            end
        end
    end
end
if isFailure == 0
    disp('FastPros test succeeded.');
else
    disp('FastPros test failed.');
end
