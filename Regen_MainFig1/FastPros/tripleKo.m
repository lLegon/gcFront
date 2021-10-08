function [TripleKoSolution, TripleKoStatus, TripleKoResult] = ...
        tripleKo(model, biomassRxn, targetRxn, oxygenRxn, options)
%tripleKo  is a function of triple knockout simulation which
%attempts to identify sets of knockout candidates for bioproductions.
% 
% [TripleKoSolution, TripleKoStatus, TripleKoResult] = ...
%         tripleKo(model, biomassRxn, targetRxn, oxygenRxn, options)
% 
% INPUTS
%  model        Structure containing all necessary variables to describe
%               a stoishiometricc model
%   rxns                    Rxns in the model
%   mets                    Metabolites in the model
%   S                       Stoichiometric matrix (sparse)
%   b                       RHS of Sv = b (usually zeros)
%   c                       Objective coefficients
%   lb                      Lower bounds for fluxes
%   ub                      Upper bounds for fluxes
%   rev                     Reversibility of fluxes
%   rxnAssociations         Each row represents reactions of original model 
%                           used for the reaction of the reduced model
%   rxnAssocMat             Matrix of reaction associations between reduced and original model
%                           (row: reactions in reduced model, column: reactions in original model)
%   unqGeneSetsMat          Matrix of genes-geneSets associations
%   geneSets                List of gene sets created by reaction combination
%   geneSetRxnMat           Matrix of geneSets-rxns associations
%                           (row: geneSets, column, rxns)
%   geneSetAssocRxns        List of reactions associated with gene sets
%   essentialGeneSets       Essential gene sets for the cell growth or the target producton
% 
%  biomassRxn       Reaction representing biomass objective function
%  targetRxn        Reaction whose flux is to be maximized
%  oxygenRxn        Reaction representing oxygen uptake
% 
% 
% OPTIONAL INPUTS
%  options
%   rxnList             Reaction list as knockout candidates 
%                       (default: all reactions in the model)
%   maxKoNum            Maximum knockout number (default: 10)
%   verbFlag            Verbose flag (default: false)
% 
% 
% OUTPUTS
%  TripleKoSolution     Solution structure
%   flux                        Each column represents flux distributions in the corresponding knockout strains
%   koGeneSetIDs                Each row represents IDs of knocked out gene sets
%   koGeneSetNum                Knockout number of gene sets
%   koGeneSets                  Each row represents knocked out gene sets
%   koRxnSets                   Each row represents reactions of knocked out gene sets
%   prodRates                   Each row represents production rates of the target compound 
%                               in the corresponding knockout strains
% 
%  TripleKoStatus       Status structure
%   existEssentialKoStrains     Column i in represents whether essential knockout strains exist 
%                               in the i th generation
%   existProdKoStrains          Column i in represents whether knockout strains with target production exist 
%                               in the i th iteration
%   possibleKoGeneSetIDs        Gene set IDs of knockout candidates
%   possibleKoRxnSets           Reaction sets encoded by gene sets of knockout candidates
%   maxKoNum                    Maximum knockout number
%   maxProdStrain               Structure of the identified strain with the maximum target production rate
%    flux                        Flux distribution of the identified strain with maximum target production rate
%    koGeneSetIDs                IDs of knocked out gene sets
%    koGeneSetNum                Knockout number of gene sets
%    koGeneSets                  Knocked out gene Sets
%    koRxnSets                   Reactions of knocked out gene sets
%    prodRate                    Production rate of the target compound
%   model                       Cobra model
%   selStrainNum                Number of strains selected as parent strains for the next iteration
%   targetInfo                  Information of target compound and reaction
%    met                         Target compound abbreviation
%    metID                       Target compound ID
%    metMW                       Molecular weight of the target compound
%    metName                     Name of target compound
%    rxn                         Target reaction whose flux is to be maximized
%    rxnID                       Target reaction ID
%   targetRxn                   Target reaction whose flux is to be maximized
%   targetRxnID                 Target reaction ID
%   theorFlux                   Theoretical flux distribution for maximized target production
%   theorProdRate               Theoretical maximum rate of target production
%   time                        Column i represents cumulative CPU time until the end of the i th generation 
% 
%  FastProsResult       Result structure of FastPros
%   essentialGeneSetCmbs        Each row represents essential gene set combinations 
%                               to reach minimum cell growth threshold or target production
%   koGeneSetNum                Total knockout number in a strain
%   prodKoStrains               Same structure to solution
% 
% 4/23/2013   Satoshi OHNO


%% Default values

tic
if nargin < 4
    error('Model, biomassRxn, targetRxn and oxygenRxn must be specified.')
end

if ~exist('options','var')
    loadData = load([cd '\' model.description '.mat'],'model');
    options.rxnList = loadData.model.rxns;
    options.maxKoNum = 3;
%     options.selStrainNum = [];
%     options.selIncUtargetStrains = true;
    options.checkSamplesNum = 1000;
    options.verbFlag = false;
    clear loadData
end
if ~isfield(options,'rxnList')
    loadData = load([cd '\' model.description '.mat'],'model');
    options.rxnList = loadData.model.rxns;
    clear loadData
end
if ~isfield(options,'maxKoNum')
    options.maxKoNum = 3;
elseif options.maxKoNum <=0 || options.maxKoNum > 3
    error('options.maxKoNum must be larger than zero and be equal to or smaller than 3.')
end
% if ~isfield(options,'selStrainNum')
%     options.selStrainNum = []; 
% elseif options.selStrainNum <= 0
%     error('options.selStrainNum must be larger than zero')
% end
% if ~isfield(options,'selIncUtargetStrains'), options.selIncUtargetStrains = true; end
if ~isfield(options,'checkSamplesNum')
    options.checkSamplesNum = 1000;
elseif options.checkSamplesNum <=0
    error('options.checkSamplesNum must be larger than zero.')
end
if ~isfield(options,'verbFlag'), options.verbFlag = false; end

biomassRxnID = findRxnIDs(model,biomassRxn);
if biomassRxnID == 0, error([num2str(biomassRxn) ' not found.']);end
targetRxnID = findRxnIDs(model,targetRxn);
if targetRxnID == 0, error([num2str(targetRxn) ' not found.']);end
oxygenRxnID = findRxnIDs(model,oxygenRxn);
if oxygenRxnID == 0, error([num2str(oxygenRxn) ' not found.']);end

metMWs = computeMW(model,model.mets,0);
targetInfo.rxn = targetRxn;
targetInfo.rxnID = targetRxnID;
tempTargetMet = targetRxn(4:end);
if strcmp(model.mets{1}(end),']')
    tempTargetMet(strfind(tempTargetMet,'(')) = '[';
    tempTargetMet(strfind(tempTargetMet,')')) = ']';
elseif strcmp(model.mets{1}(end),')')
    tempTargetMet(strfind(tempTargetMet,'[')) = ')';
    tempTargetMet(strfind(tempTargetMet,'[')) = ')';
end
targetInfo.metID = find(strcmp(model.mets, tempTargetMet));
if isempty(targetInfo.metID)
    tempTargetMet(strfind(tempTargetMet,'_')) = '-';
    targetInfo.metID = find(strcmp(model.mets, tempTargetMet));
end
targetInfo.met = tempTargetMet;
targetInfo.metID = find(strcmp(model.mets, targetInfo.met));
targetInfo.metName = model.metNames{targetInfo.metID};
targetInfo.metMW = metMWs(targetInfo.metID);
targetInfo = orderfields(targetInfo);


%% FBA on no deletion strain

model = changeObjective(model,{biomassRxn,targetRxn},[1,10^-5]);
solWT = optimizeCbModel(model,'max');
solWT.x = round(solWT.x*10^6)/10^6;
if solWT.x(biomassRxnID) <= model.lb(biomassRxnID);
    error('The growth rate of wild type strain is less than growth threshold')
% else
%     if solWT.x(targetRxnID) > 0
%         error('Wild type strain produces target compounds, so FastPros cannont be applied')
%     end
end


%% Maximize target production

modelMaxTarget = changeObjective(model,{biomassRxn,targetRxn},[0,1]);
modelMaxTarget = changeRxnBounds(modelMaxTarget,{oxygenRxn},-1000,'l');
solMaxTarget = optimizeCbModel(modelMaxTarget,'max');
solMaxTarget.x = round(solMaxTarget.x*10^6)/10^6;
theorFlux = solMaxTarget.x;
theorProdRate = solMaxTarget.x(targetRxnID);
if theorProdRate == 0
    error('Target compound cannot be produced.')
end


%% Remove essential genes from knockout candidates

possibleKoGeneSetIDs = find(model.essentialGeneSets==0);

inputRxnIDs = findCombinedRxnIDs(model,options.rxnList);
inputRxnIDs = inputRxnIDs(inputRxnIDs~=0);
inputGeneSets = findGeneSetsFromRxns(model,model.rxns(inputRxnIDs));
inputGeneSetIDs = findGeneSetIDs(model,inputGeneSets);
possibleKoGeneSetIDs = intersect(possibleKoGeneSetIDs,inputGeneSetIDs);

% Convert gene set IDs to reactions
possibleKoRxnSets = findRxnSetsFromGeneSetIDs(model,possibleKoGeneSetIDs);

if size(possibleKoGeneSetIDs,1) < size(possibleKoGeneSetIDs,2)
    possibleKoGeneSetIDs = possibleKoGeneSetIDs';
    possibleKoRxnSets = possibleKoRxnSets';
end
possibleKoGeneNum = length(possibleKoGeneSetIDs);


%% Initialize TKO results

if options.verbFlag == true
    disp('Triple knockout simulation start.')
    disp(['Total ' num2str(possibleKoGeneNum) ' gene sets can be knocked out.'])
end
TripleKoStatus.maxKoNum = 0;
TripleKoStatus.maxProdStrain.koGeneSetNum = 0;
TripleKoStatus.maxProdStrain.koGeneSetIDs = 0;
TripleKoStatus.maxProdStrain.koGeneSets = [];
TripleKoStatus.maxProdStrain.koRxnSets = [];
% TKOStatus.maxProdStrain.utarget = 0;
TripleKoStatus.maxProdStrain.prodRate = 0;
TripleKoStatus.existProdKoStrains = false(1,options.maxKoNum);
TripleKoStatus.existEssentialKoStrains = false(1,options.maxKoNum);
existProdKoStrains = false(1,options.maxKoNum);
existEssentialKoStrains = false(1,options.maxKoNum);
time = zeros(1,options.maxKoNum+1);
startTime = toc;


%% Simulate wild type strain

flux0 = glpk(-model.c,model.S,model.b,model.lb,model.ub,repmat('S',length(model.mets),1));
time(1) = toc-startTime;

possibleKoGeneNum = length(possibleKoGeneSetIDs);

isKnockedOut=zeros(possibleKoGeneNum,1);
for i = 1 : possibleKoGeneNum
    tempRxnTFs = logical(model.geneSetRxnMat(possibleKoGeneSetIDs(i),:));
    isKnockedOut(i) = sum(flux0(tempRxnTFs) ~= 0,1);
end
zeroFluxGeneSetIDs0 = possibleKoGeneSetIDs(isKnockedOut==0);  % 野生株で0のフラックスを持つ遺伝子


%% Simulate knockout strains

fbaLpNum = 1;
for tempKoNum = 1:options.maxKoNum
    if options.verbFlag == true
        disp(' ')
        disp([num2str(tempKoNum) ' KO'])
    end

    switch tempKoNum
        case 1
            koGeneSetIDs = possibleKoGeneSetIDs;
        case 2
            koGeneSetIDs = zeros(nchoosek(length(possibleKoGeneSetIDs),2),2);
            ii = 0;
            for i1 = 2: length(possibleKoGeneSetIDs)
                koGeneSetIDs(ii+1:ii+i1-1,1) = possibleKoGeneSetIDs(i1);
                koGeneSetIDs(ii+1:ii+i1-1,2) = possibleKoGeneSetIDs(1:i1-1);
                ii = ii + i1-1;
            end
        case 3
            koGeneSetIDs = zeros(nchoosek(length(possibleKoGeneSetIDs),3),3);
            ii = 0;
            for i2 = 3: length(possibleKoGeneSetIDs)
                for i1 = 2: i2-1
                    koGeneSetIDs(ii+1:ii+i1-1,1) = possibleKoGeneSetIDs(i2);
                    koGeneSetIDs(ii+1:ii+i1-1,2) = possibleKoGeneSetIDs(i1);
                    koGeneSetIDs(ii+1:ii+i1-1,3) = possibleKoGeneSetIDs(1:i1-1);
                    ii = ii + i1-1;
                end
            end
    end
    
    if tempKoNum == 1
        hasNoEffects = ismember(possibleKoGeneSetIDs,zeroFluxGeneSetIDs0);        
        growth = zeros(size(koGeneSetIDs,1),1);
        growth(hasNoEffects) = -888;
        prodRates = zeros(size(koGeneSetIDs,1),1);        
    else
        hasNoEffects = ismember(sort(koGeneSetIDs,2),sort(zeroFluxGeneSetIDs,2),'rows');        
        growth = zeros(size(koGeneSetIDs,1),1);
        growth(hasNoEffects) = -888;
        prodRates = zeros(size(koGeneSetIDs,1),1);        
    end
    clear zeroFluxGeneSetIDs
    
    simulationList = find(growth == 0);
    if options.verbFlag
        disp([num2str(length(simulationList)) ' times FBA calculation.'])
    end
    
    isKoGenes = sparse(length(simulationList),length(model.geneSets));
    for i = 1 : length(simulationList)
        isKoGenes(i,koGeneSetIDs(simulationList(i),:)) = 1;
    end
    tempIsKoRxns = isKoGenes * model.geneSetRxnMat > 0;
    
    parSol = cell(1,length(simulationList));
    parfor i = 1 : length(simulationList)
        modelDel = model;
        modelDel.ub(tempIsKoRxns(i,:),1)=0;
        modelDel.lb(tempIsKoRxns(i,:),1)=0;
        [x,f,~,~]= ...
            glpk(modelDel.c,modelDel.S,modelDel.b,modelDel.lb,modelDel.ub,...
            repmat('S',length(modelDel.mets),1),repmat('C',length(modelDel.rxns),1),-1);
        parSol{i}.x = x;
        parSol{i}.f = f;
    end
    clear tempIsKoRxns
    sol = [parSol{:}];
    clear parSol
    
    flux = sparse(round(cat(2,[sol.x])*10^6)/10^6);
    objv = cat(2,[sol.f]);
    fbaLpNum = fbaLpNum + length(simulationList);
    
    essentialIDs = find(flux(biomassRxnID,:)<=model.lb(biomassRxnID));
    flux(:,essentialIDs) = 0;
    objv(essentialIDs) = 0;
    growth(simulationList)=objv';
    prodRates(simulationList) = flux(targetRxnID,:);
    clear sol
    
    isKnockedOut=zeros(possibleKoGeneNum,length(simulationList));
    for i = 1 : possibleKoGeneNum
        tempRxnTFs = logical(model.geneSetRxnMat(possibleKoGeneSetIDs(i),:));
        isKnockedOut(i,:) = sum(flux(tempRxnTFs,:)~=0,1);
    end
    [noEffectsKoGeneSetIDs,sameStrainIDs] = find(isKnockedOut==0);
    zeroFluxGeneSetIDs = [koGeneSetIDs(simulationList(sameStrainIDs),:),possibleKoGeneSetIDs(noEffectsKoGeneSetIDs)];
    switch tempKoNum
        case 1
            zeroFluxGeneSetIDs1 = zeroFluxGeneSetIDs;
            zeroFluxGeneSetIDs = [zeroFluxGeneSetIDs;nchoosek(zeroFluxGeneSetIDs0,2)];
        case 2
            zeroFluxGeneSetIDs = [zeroFluxGeneSetIDs;nchoosek(zeroFluxGeneSetIDs0,3)];
            [unique1KoGeneSets,uniqueFirstID] = unique(zeroFluxGeneSetIDs1(:,1),'first');
            [~,uniqueLastID] = unique(zeroFluxGeneSetIDs1(:,1),'last');
            parZeroFluxGeneSetIDs = cell(length(uniqueFirstID),1);
            for i = 1 : length(uniqueFirstID)
                temp = zeros(nchoosek(uniqueLastID(i)-uniqueFirstID(i)+1,2),3);
                temp(:,1) = unique1KoGeneSets(i);
                temp(:,2:3) = nchoosek(zeroFluxGeneSetIDs1(uniqueFirstID(i):uniqueLastID(i),2),2);
                parZeroFluxGeneSetIDs{i} = temp;
            end
            zeroFluxGeneSetIDs = [zeroFluxGeneSetIDs;cat(1,parZeroFluxGeneSetIDs{:})];
            clear parZeroFluxGeneSetIDs            
    end
%     zeroFluxGeneSetIDs = unique(sort(zeroFluxGeneSetIDs,2),'rows');

    
    %Integrate essential gene set combinations to reach minimun cell growth
%     essentialIDs = find(flux(biomassRxnID,:)<=model.lb(biomassRxnID));
    if isempty(essentialIDs)
        existEssentialKoStrains(tempKoNum) = true;
        essentialGeneSetCombs = [];
    else
        essentialGeneSetCombs = koGeneSetIDs(simulationList(essentialIDs),:);
    end
    
    %Identify target produciton strains
    prodIDs = find(flux(targetRxnID,:) > 0);
    if isempty(prodIDs)
        prodKoStrains = struct([]);
        if options.verbFlag == true
            disp('No production strains were identified.')
        end
    else
        prodKoStrains.koNum = tempKoNum;
        prodKoStrains.koGeneSetIDs = koGeneSetIDs(simulationList(prodIDs),:);
        prodKoStrains.koRxnSets = findRxnSetsFromGeneSetIDs(model,prodKoStrains.koGeneSetIDs);
        prodKoStrains.prodRates = full(flux(targetRxnID,prodIDs)');
        prodKoStrains.flux = full(flux(:,prodIDs));
        
        [~,sortID] = sort(prodKoStrains.prodRates,'descend');
        prodKoStrains.koGeneSetIDs = prodKoStrains.koGeneSetIDs(sortID,:);
        prodKoStrains.koGeneSets = model.geneSets(prodKoStrains.koGeneSetIDs)';
        prodKoStrains.koRxnSets = prodKoStrains.koRxnSets(sortID,:);
        prodKoStrains.prodRates = prodKoStrains.prodRates(sortID);
        prodKoStrains.flux = prodKoStrains.flux(:,sortID);
        
        existProdKoStrains(tempKoNum) = true;
        
        maxID = 1;
        tempMaxProdRate = prodKoStrains.prodRates(maxID);
        if options.verbFlag == true
            disp([num2str(length(prodKoStrains.prodRates)) ' production strains were identified.'])
            disp(['Maximum target production rate is ' num2str(tempMaxProdRate) ' mmol/gDW/h '])
%             disp(['utarget = ' num2str(productionKoStrains.utarget(maxID,tempKoNum))])
            for i = 1:tempKoNum
                disp([''  prodKoStrains.koGeneSets{maxID,i} ' gene set is knocked out.'])
                disp([' ' prodKoStrains.koRxnSets{maxID,i} ' reaction set is knocked out.'])
            end
        end
    end

    
    tempResult.essentialGeneSetCombinations = essentialGeneSetCombs;
    tempResult.koNum = tempKoNum;
    tempResult.prodKoStrains = prodKoStrains;
    TripleKoResult(tempKoNum) = tempResult;
    clear productionKoStrains
    time(tempKoNum+1) = toc-startTime;
    
    if options.verbFlag == true
        toc
    end
    clear prodKoStrains
end


%% Arrange results

proNum = 0;
if any(existProdKoStrains)
    for i = 1:options.maxKoNum
        if ~isempty(TripleKoResult(i).prodKoStrains)
            proNum = proNum + 1;
            TripleKoSolution(proNum) = TripleKoResult(i).prodKoStrains;
            TripleKoSolution(proNum).flux = full(TripleKoSolution(proNum).flux);
        end
    end
%     TripleKoSolution.flux = full(TripleKoSolution.flux);
    tempMaxProdRate = zeros(1,length(TripleKoSolution));
    for i = 1 : length(TripleKoSolution)
        tempMaxProdRate(i) = TripleKoSolution(i).prodRates(1);
    end
    [maxProdRate,maxID] = max(tempMaxProdRate);
    maxProdKoNum = TripleKoSolution(maxID).koNum;

    maxProdStrain.koNum =  TripleKoSolution(maxID).koNum;
    maxProdStrain.koGeneSetIDs = TripleKoSolution(maxID).koGeneSetIDs(1,:);
    maxProdStrain.koGeneSets = TripleKoSolution(maxID).koGeneSets(1,:);
    maxProdStrain.koRxnSets = TripleKoSolution(maxID).koRxnSets(1,:);
    maxProdStrain.prodRate = maxProdRate;
    maxProdStrain.flux = TripleKoSolution(maxID).flux(:,1);

else
    TripleKoSolution = struct([]);
    maxProdRate = 0;
    maxProdKoNum = 0;
end

TripleKoStatus.time = time;
TripleKoStatus.existProdKoStrains = existProdKoStrains;
TripleKoStatus.maxKoNum = options.maxKoNum;
TripleKoStatus.maxProdKoNum = maxProdKoNum;
TripleKoStatus.maxProdRate = maxProdRate;
TripleKoStatus.fbaLpNum = fbaLpNum;
TripleKoStatus.possibleKoGeneSetIDs = possibleKoGeneSetIDs;
TripleKoStatus.possibleKoRxnSets = possibleKoRxnSets;
TripleKoStatus.maxProdStrain = maxProdStrain;
TripleKoStatus.model = model;
TripleKoStatus.theorFlux = theorFlux;
TripleKoStatus.theorProdRate = theorProdRate;
TripleKoStatus.targetRxn = targetRxn;
TripleKoStatus.targetRxnID = targetRxnID;
TripleKoStatus.targetInfo = targetInfo;
TripleKoStatus.existProdKoStrains = existProdKoStrains;
TripleKoStatus.existEssentialKoStrains = existEssentialKoStrains;

TripleKoSolution = orderfields(TripleKoSolution);
TripleKoStatus = orderfields(TripleKoStatus);


end


