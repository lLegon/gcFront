function [modelReduced,biomassRxn,targetRxn,oxygenRxn]= ...
    reduceModelForFP(model,biomassRxn,targetRxn,oxygenRxn,options)
%reduceModelForFP  is a function to reduce a COBRA metabolic model for FastPros
%
% [modelNew,biomassRxn,targetRxn,oxygenRxn]= ...
%     reduceModelForFP(model,biomassRxn,targetRxn,oxygenRxn,options)
% 
% INPUTS
% model         Structure containing following required fields to describe a COBRA stoichiometric model
%   rxns                    Reaction name abbreviation; reaction ID; order corresponds to S matrix.
%   mets                    Metabolite name abbreviation; metaboliteID; order corresponds to S matrix
%   S                       Stoichiometric matrix in sparse format
%   b                       RHS of Sv = b (usually zeros)
%   c                       Objective coefficient for corresponding reactions
%   lb                      Lower flux bound for corresponding reactions
%   ub                      Upper flux bound for corresponding reactions
%   rev                     Logical array; true for reversible reactions, otherwise false
%   genes                   List of all genes within the model
%   rxnGeneMat              Matrix with rows corresponding to reactions and columns corresponding to genes
%   grRules                 Rules field in a format readable format
%   metFormulas             Elemental formula for each metabolite
%                 
%  biomassRxn               Reaction reresenting biomass objective function
%  targetRxn                Reaction whose flux is to be maximized
%  oxygenRxn                Reaction representing oxygen uptake
% 
% OPTIONAL INPUTS
%  options
%   verbFlag                Verbose flag
%   loadFVAFlux             Load maximum and minimun fluxes of each reaction calculated by flux variability analysis
%                           (default: false)
% 
% OUTPUTS
%  modelReduced      COBRA model structure added with the following generated fields.
%                    In the filed "rxns", Combined reactions are represented using "/", as "reactionA/reactionB". 
%    rxnFormulas            Reaction formulas
%    rxnAssociations        Reactions drived from original model
%    rxnAssocMat            Matrix of reaction associations between reduced and original model
%                           (row: rxns in reduced model, column: rxns in original model)
%    unqGeneSetsMat         Matrix of genes-geneSets associations
%    geneSets               List of gene sets
%    geneSetRxnMat          Matrix of geneSets-rxns associations
%                           (row: geneSets, column, rxns)
%    geneSetAssocRxns       List of reactions associated with gene sets
%    essentialGeneSets      Essential gene sets for the cell growth or the target producton
%    oriRxns                Reactions in the original model
%    reductionStatus        Status representing whether model reduction was success or not.
%                           1: Model reduction was success.
%                           2: Growth rate of wild type strain was changed by model reduction.  
%                           3: Growth rate of single knockout strains were changed by model reduction.
%   
%  biomassRxn       Reaction representing biomass objective function
%  targetRxn        Reaction whose flux is to be maximized
%  oxygenRxn        Reaction representing oxygen uptake
% 
% 
% Aug. 5th, 2013    Satoshi OHNO


tic

%% Check if the model was already reduced or not

if isfield(model, 'reductionStatus')
    modelReduced = model;
    warning('The model was already reduced')
    return
end

%% Prepare variables

if ~exist('options','var') || isempty(options)
    options.verbFlag = false;
    options.loadFVAFlux = false;
else
    if ~isfield(options,'verbFlag'), options.verbFlag = false; end
    if ~isfield(options,'loadFVAFlux'), options.loadFVAFlux = false; end
end

[~,n] = size(model.S);
biomassRxnID = findRxnIDs(model,biomassRxn);
if biomassRxnID == 0, error([num2str(biomassRxn) ' not found.']);end
targetRxnID = findRxnIDs(model,targetRxn);
if targetRxnID == 0, error([num2str(targetRxn) ' not found.']);end
oxygenRxnID = findRxnIDs(model,oxygenRxn);
if oxygenRxnID == 0, error([num2str(oxygenRxn) ' not found.']);end

rxnFormulas = printRxnFormula(model, model.rxns, false);
model.rxnFormulas = rxnFormulas;
reductionStatus = 1;


%% FBA on no deletion strain

model = changeObjective(model,{biomassRxn,targetRxn},[1,10^-5]);
solWTx = solveModel(model);
solWTx(isnan(solWTx))=0;
solWTx = round(solWTx*10^6)/10^6;
if solWTx(biomassRxnID) <= model.lb(biomassRxnID);
    error('The growth rate of wild type strain is less than growth threshold')
else
    if solWTx(targetRxnID) > 0
        error('Wild type strain produces target metabolites, so FastPros cannont be applied')
    end
end


%% Calculate theoretical maximum production rate of target metabolite

modelMaxTarget = changeObjective(model,{biomassRxn,targetRxn},[0,1]);
modelMaxTarget = changeRxnBounds(modelMaxTarget,{oxygenRxn},-1000,'l');
solMaxTargetx = solveModel(modelMaxTarget);
solMaxTargetx(isnan(solMaxTargetx))=0;
solMaxTargetx = round(solMaxTargetx*10^6)/10^6;
theorMaxProdRate = solMaxTargetx(targetRxnID);
if theorMaxProdRate == 0
    error('Target metabolite cannot be produced.')
end

%% Flux Variability analysis

if options.verbFlag == true
    disp('Model reduction start.')
end

if options.loadFVAFlux == true
    load(['fluxFVA_' model.description '.mat']);
else
    if options.verbFlag == true
        disp('Flux variability anaysis')
    end
    parFluxFVA = cell(length(model.rxns),1);
    for i = 1 : length(model.rxns)
        modelFVA = model;
        modelFVA.c = zeros(length(model.rxns),1);
        modelFVA.c(i) = 1;
%         solMin = glpk(modelFVA.c,modelFVA.S,modelFVA.b,modelFVA.lb,modelFVA.ub,repmat('S',length(modelFVA.mets),1));
        solMin = solveModel(modelFVA);
        modelFVA.c = zeros(length(model.rxns),1);
        modelFVA.c(i) = -1;
%         solMax = glpk(modelFVA.c,modelFVA.S,modelFVA.b,modelFVA.lb,modelFVA.ub,repmat('S',length(modelFVA.mets),1));
        solMax = solveModel(modelFVA);
        parFluxFVA{i} = [solMin(i),solMax(i)];
    end
    fluxFVA = cat(1,parFluxFVA{:});
    fluxFVA(isnan(fluxFVA))=0;
%     fluxFVA = zeros(length(model.rxns),2);
%     [fluxFVA(:,1),fluxFVA(:,2)]= fluxVariability(model,0);
    fluxFVA = round(fluxFVA*10^6)/10^6;
%     save(['fluxFVA_' model.description '.mat'],'fluxFVA')
end


%% Remove reactions which cannot carry any fluxes in the given condition

if options.verbFlag == true
    disp('Remove reactions from model')
end
removeRxnIDs = find(sum(abs(fluxFVA),2)==0);
modelRem = removeRxns(model,model.rxns(removeRxnIDs));

% Move biomassRxn, targetRxn, and oxygenRxn to first
minLb = min(model.lb);
maxUb = max(model.ub);
specifiedRxnIDs = find(...
    ismember(model.lb,[0,minLb])==0 | ismember(model.ub,[0,maxUb])==0);
% specifiedRxnIDs = [biomassRxnID;targetRxnID;oxygenRxnID;setdiff(specifiedRxnIDs,[biomassRxnID,targetRxnID,oxygenRxnID])'];
specifiedRxnIDs = [biomassRxnID;targetRxnID;oxygenRxnID;setdiff(specifiedRxnIDs,[biomassRxnID,targetRxnID,oxygenRxnID])];
specifiedRxnIDs = [specifiedRxnIDs(1:3);setdiff(specifiedRxnIDs(4:end),removeRxnIDs)];
specifiedRxnLbs = model.lb(specifiedRxnIDs);
specifiedRxnUbs = model.ub(specifiedRxnIDs);
specifiedRxnsMat = zeros(n,length(specifiedRxnIDs));
for i = 1 : length(specifiedRxnIDs)
    specifiedRxnsMat(specifiedRxnIDs(i),i) = 1;
end
specifiedRxnsMat = specifiedRxnsMat(setdiff(1:n,removeRxnIDs),:);
clear rxnRemoveNum

nonRemovedGene = sum(modelRem.rxnGeneMat,1)>0;
modelRem.rxnGeneMat = modelRem.rxnGeneMat(:,nonRemovedGene);
modelRem.genes = modelRem.genes(nonRemovedGene);
clear nonRemovedGene

if options.verbFlag == true
    disp([num2str(length(modelRem.rxns)) ' reactions' ...
        ', ' num2str(length(modelRem.mets)) ' metabolites, ' num2str(length(modelRem.genes)) ' genes'])
end


%% Combine adjacent reactions without branching

if options.verbFlag == true
    disp('Combine adjacent reactions without branching')
end

% Identify exchange reactions and its flux bounds
[isExRxns, isUptRxns] = findExcRxns(modelRem);

% Combine linear reactions (except reactions with metabolites in biomassRxn)
combineRest = 1;
rxnSets = modelRem.rxns;
tempBiomassRxnID = findRxnIDs(modelRem,biomassRxn);
tempTargetRxnID = findRxnIDs(modelRem,targetRxn);
nonBiomassTargetMetIDs = find(modelRem.S(:,tempBiomassRxnID)==0 & modelRem.S(:,tempTargetRxnID)==0); 
while combineRest == 1
    combineRest = 0;
    
    for tempMetID = nonBiomassTargetMetIDs'
        metAssocRxnIDs = find(modelRem.S(tempMetID,:));
        
        if size(metAssocRxnIDs ,2) == 2
            combine = 0;
            
            switch sum(modelRem.rev(metAssocRxnIDs))
                case 2   %Both reactions are reversible 
                    combineRest = 1;
                    combine = 1;
                    if isExRxns(metAssocRxnIDs(2)) >= 1
                        temp = metAssocRxnIDs;
                        metAssocRxnIDs(1) = temp(2);
                        metAssocRxnIDs(2) = temp(1);                        
                    end
                case 1  %One of reactions is reversible
                    combineRest = 1;
                    combine = 1;
                    if modelRem.rev(metAssocRxnIDs(1)) == 1
                        temp = metAssocRxnIDs;
                        metAssocRxnIDs(1) = temp(2);
                        metAssocRxnIDs(2) = temp(1);
                        clear temp
                    end
                case 0  %Both reactions are irreversible
                    if prod(modelRem.S(tempMetID , metAssocRxnIDs)) < 0
                        combineRest = 1;
                        combine = 1;
                        if isExRxns(metAssocRxnIDs(2)) >= 1
                            temp = metAssocRxnIDs;
                            metAssocRxnIDs(1) = temp(2);
                            metAssocRxnIDs(2) = temp(1);
                        end
                    end
            end
            
            if combine == 1
                
                %Convert stoichiometric matrix
                coeffRatio = modelRem.S(tempMetID,metAssocRxnIDs(1))/modelRem.S(tempMetID,metAssocRxnIDs(2));
                modelRem.S(:,metAssocRxnIDs(2)) = modelRem.S(:,metAssocRxnIDs(2)) * (-coeffRatio);
                modelRem.S(:,metAssocRxnIDs(1)) = ...
                    modelRem.S(:,metAssocRxnIDs(1)) + modelRem.S(:,metAssocRxnIDs(2));
                modelRem.S(:,metAssocRxnIDs(2)) = 0 ;
                
                %Convert objective function
                modelRem.c(metAssocRxnIDs(2)) = modelRem.c(metAssocRxnIDs(2)) * (-coeffRatio);
                modelRem.c(metAssocRxnIDs(1)) = ...
                    modelRem.c(metAssocRxnIDs(1)) + modelRem.c(metAssocRxnIDs(2));
                modelRem.c(metAssocRxnIDs(2)) = 0 ;            
                
                %Convert reaction and metabolite names
                modelRem.rxns(metAssocRxnIDs(1)) = ...
                    {[modelRem.rxns{metAssocRxnIDs(1)} , '/' , modelRem.rxns{metAssocRxnIDs(2)}]};
                if isfield(modelRem,'rxnNames')
                    modelRem.rxnNames(metAssocRxnIDs(1)) = ...
                        {[modelRem.rxnNames{metAssocRxnIDs(1)} , '/' , modelRem.rxnNames{metAssocRxnIDs(2)}]};
                    modelRem.rxnNames(metAssocRxnIDs(2)) = {'none'};
                end
                rxnSets(metAssocRxnIDs(1)) = ...
                    {[rxnSets{metAssocRxnIDs(1)} , '/' , rxnSets{metAssocRxnIDs(2)}]};
                rxnSets(metAssocRxnIDs(2)) = {'none'};
                modelRem.mets(tempMetID,1) = {'none'};
                
                %Convert genes
                modelRem.rxnGeneMat(metAssocRxnIDs(1),:) = ...
                    modelRem.rxnGeneMat(metAssocRxnIDs(1),:) + modelRem.rxnGeneMat(metAssocRxnIDs(2),:);
                modelRem.rxnGeneMat(metAssocRxnIDs(2),:) = 0;
                modelRem.grRules(metAssocRxnIDs(1),1) = ...
                    {[modelRem.grRules{metAssocRxnIDs(1)} , '/' , modelRem.grRules{metAssocRxnIDs(2)}]};
                modelRem.grRules(metAssocRxnIDs(2),1) = {'none'};
                                
                %Convert exchange reactions
                if isExRxns(metAssocRxnIDs(2)) > 0
                    isExRxns(metAssocRxnIDs(1)) =isExRxns(metAssocRxnIDs(1)) + isExRxns(metAssocRxnIDs(2)) ;
                    isExRxns(metAssocRxnIDs(2)) = 0;
                    if isUptRxns(metAssocRxnIDs(2)) > 0
                        isUptRxns(metAssocRxnIDs(1)) = ...
                            isUptRxns(metAssocRxnIDs(1)) + isUptRxns(metAssocRxnIDs(2)) ;
                        isUptRxns(metAssocRxnIDs(2)) = 0;
                    end
                end
                
                %Convert specified reactions
                specifiedRxnsMat(metAssocRxnIDs(1),:) = ...
                    specifiedRxnsMat(metAssocRxnIDs(1),:) + specifiedRxnsMat(metAssocRxnIDs(2),:) * (-coeffRatio);
                specifiedRxnsMat(metAssocRxnIDs(2),:) = 0;
                                
            end
        end 
        clear metAssocRxnsNum coefRatio        
    end
end
clear combine tempBiomassRxnNum nonBiomassMets coeffRatio

combineRxnIDs = find((sum(modelRem.S.^2) == 0));
modelCom = removeRxns(modelRem,modelRem.rxns(combineRxnIDs));
nonCombineRxnIDs = setdiff(1:length(modelRem.rxns),combineRxnIDs); 
isExRxns = isExRxns(nonCombineRxnIDs);
isUptRxns = isUptRxns(nonCombineRxnIDs);
specifiedRxnsMat = specifiedRxnsMat(nonCombineRxnIDs,:);
clear rxnCombineNum

modelCom.rxnGeneMat(modelCom.rxnGeneMat>=1) = 1;
nonCombinedGene = find(sum(modelCom.rxnGeneMat,1)>0);
modelCom.rxnGeneMat = modelCom.rxnGeneMat(:,nonCombinedGene);
modelCom.genes = modelCom.genes(nonCombinedGene);
clear nonCombinedGene

if options.verbFlag == true
    disp([num2str(length(modelCom.rxns)) ' reactions' ...
        ', ' num2str(length(modelCom.mets)) ' metabolites, ' num2str(length(modelCom.genes)) ' genes'])
end


%% Move exchange reactions to forward

[specifiedRxnIDs,~] = find(specifiedRxnsMat);
moveRxnIDs = columnVector([specifiedRxnIDs; setdiff(find(isExRxns),specifiedRxnIDs)]);
modelMove = moveRxnFirst(modelCom,moveRxnIDs);
nonMoveRxnIDs = columnVector(setdiff(1:length(modelCom.rxns),moveRxnIDs));
oldExflag = isExRxns;
isExRxns = oldExflag([moveRxnIDs;nonMoveRxnIDs]);
oldMinMedium = isUptRxns;
isUptRxns = oldMinMedium([moveRxnIDs;nonMoveRxnIDs]);
oldSpecifiedRxnsMat = specifiedRxnsMat;
specifiedRxnsMat = oldSpecifiedRxnsMat([moveRxnIDs;nonMoveRxnIDs],:);


%% Convert direction of exchange reactions

modelMoveMetMWs = computeMW(modelMove,modelMove.mets,0);
for i = length(find(specifiedRxnsMat))+1: length(find(specifiedRxnsMat))+length(find(isExRxns))
    leftMetsNum = find(modelMove.S(:,i) < 0);
    rightMetsNum = find(modelMove.S(:,i) > 0);
    if isempty(leftMetsNum); leftMWs=0;
    else leftMWs = -sum(modelMoveMetMWs(leftMetsNum) .* modelMove.S(leftMetsNum,i));
    end
    if isempty(rightMetsNum); rightMWs=0;
    else rightMWs = sum(modelMoveMetMWs(rightMetsNum) .* modelMove.S(rightMetsNum,i));
    end
    
    %Conver direction of reversible exchnage reactions if sum of molecular
    %weight in right hands is larger than that in left hands
    if rightMWs - leftMWs >= 2 
        if modelMove.rev(i) == 0 
            modelMove.S(:,i) = - modelMove.S(:,i);
            specifiedRxnsMat(i,:) = -specifiedRxnsMat(i,:);
        end
    end
end


%% Create reaction formulas of reduced model

rxnFormulas = printRxnFormula(modelMove, modelMove.rxns, false);


%% Identify reaction correlations of models between before and after reduction

rxnAssociations = cell(length(modelMove.rxns),100);
rxnAssocMat = zeros(length(modelMove.rxns),length(model.rxns));
combineNum = zeros(length(modelMove.rxns),1);
for i = 1 : length(modelMove.rxns)
        slashLocations = find(modelMove.rxns{i,1} == '/');
        combineNum(i) = length(slashLocations)+1;
        switch length(slashLocations)
            case 0
                rxnAssociations(i,1)=modelMove.rxns(i,1);
                rxnAssocMat(i,strcmp(modelMove.rxns(i),model.rxns))= 1;
            case 1
                rxnAssociations(i,1) = {modelMove.rxns{i}(1:slashLocations(1)-1)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(1:slashLocations(1)-1),model.rxns))= 1;
                rxnAssociations(i,2) = {modelMove.rxns{i}(slashLocations(1)+1:end)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(slashLocations(1)+1:end),model.rxns))= 1;
            otherwise
                rxnAssociations(i,1) = {modelMove.rxns{i}(1:slashLocations(1)-1)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(1:slashLocations(1)-1),model.rxns))= 1;
                for j = 2 : length(slashLocations)
                    rxnAssociations(i,j) = {modelMove.rxns{i}(slashLocations(j-1)+1:slashLocations(j)-1)};
                    rxnAssocMat(i,strcmp(modelMove.rxns{i}(slashLocations(j-1)+1:slashLocations(j)-1),model.rxns))= 1;
                end
                rxnAssociations(i,j+1) = {modelMove.rxns{i}(slashLocations(j)+1:end)};
                rxnAssocMat(i,strcmp(modelMove.rxns{i}(slashLocations(j)+1:end),model.rxns))= 1;
        end
end
rxnAssociations = rxnAssociations(:,1:max(combineNum));
rxnAssocMat = sparse(logical(rxnAssocMat));
clear slashLocations


%% Identify correlations between gene sets and reactions in reduced model

unqGeneSetsMat = unique(modelMove.rxnGeneMat,'rows');
geneSets = cell(size(unqGeneSetsMat,1),1);
for i = 1 : size(unqGeneSetsMat,1)
    geneNum = find(unqGeneSetsMat(i,:));
    switch length(geneNum)
        case 0
            geneSets{i} = 'None or Unknown';
        case 1
            geneSets{i} = modelMove.genes{geneNum};
        case 2
            geneSets{i} = [modelMove.genes{geneNum(1)} '_' modelMove.genes{geneNum(2)}];
        otherwise
            tempGeneSets = modelMove.genes{geneNum(1)};
            for j = 1 : length(geneNum)-1
                tempGeneSets = [tempGeneSets '_' modelMove.genes{geneNum(j+1)}];
            end
            geneSets{i} = tempGeneSets;
    end
end

[~, rxnGeneSets] = ismember(modelMove.rxnGeneMat,unqGeneSetsMat,'rows');
geneSetRxnMat = zeros(length(geneSets),length(modelMove.rxns));
for i = 1 : length(modelMove.rxns)
    geneSetRxnMat(rxnGeneSets(i),i) = 1;
end
geneSetRxnMat = sparse(geneSetRxnMat);
clear tempGeneSets


%% Identify reactions associated with gene sets

geneSetAssocRxns = cell(length(geneSets),1);
for i = 1 : length(geneSets)
    associatedRxns = modelMove.rxns(logical(geneSetRxnMat(i,:)));
    tempAssocRxnsName = associatedRxns{1};
    if length(associatedRxns) >= 2
        for j = 2 : length(associatedRxns)
            tempAssocRxnsName(end+1:end+length(associatedRxns{j})+1) =[',' associatedRxns{j}];
        end
    end
    geneSetAssocRxns(i) = {tempAssocRxnsName};
end


%%  Set bounds of reaction fluxes

modelMove = changeRxnBounds(modelMove,modelMove.rxns(modelMove.rev==0),0,'l');
modelMove = changeRxnBounds(modelMove,modelMove.rxns(modelMove.rev==0),maxUb,'u');
modelMove = changeRxnBounds(modelMove,modelMove.rxns(modelMove.rev==1),minLb,'l');
modelMove = changeRxnBounds(modelMove,modelMove.rxns(modelMove.rev==1),maxUb,'u');
modelMove = changeRxnBounds(modelMove,modelMove.rxns(isExRxns==1),0,'l');
modelMove = changeRxnBounds(modelMove,modelMove.rxns(isUptRxns==1),minLb,'l');

[specifiedRxnIDs,cIDs] = find(specifiedRxnsMat);
for i = 1 : length(specifiedRxnIDs)
    coef = specifiedRxnsMat(specifiedRxnIDs(i),cIDs(i));
    if coef > 0
        if specifiedRxnLbs(i) == minLb
            modelMove.lb(specifiedRxnIDs(i)) = minLb;
        else modelMove.lb(specifiedRxnIDs(i)) = specifiedRxnLbs(i) / coef;
        end
        if specifiedRxnUbs(i) == maxUb
            modelMove.ub(specifiedRxnIDs(i)) = maxUb;
        else modelMove.ub(specifiedRxnIDs(i)) = specifiedRxnUbs(i) / coef;
        end
    else
        if specifiedRxnUbs(i) == maxUb
            modelMove.lb(specifiedRxnIDs(i)) = minLb;
        else modelMove.lb(specifiedRxnIDs(i)) = specifiedRxnUbs(i) / coef;
        end
        if specifiedRxnLbs(i) == minLb
            modelMove.ub(specifiedRxnIDs(i)) = maxUb;
        else modelMove.ub(specifiedRxnIDs(i)) = specifiedRxnLbs(i) / coef;
        end
    end
end


%% Confirm model reduction by FBA simulation

% FBA simulation on wild type strain
testSolx = solveModel(modelMove); 
testSolx(isnan(testSolx))=0;
testSolx = round(testSolx*10^6)/10^6;
if options.verbFlag == true
    disp(['Biomass production rate     before model reductionF' num2str(solWTx(biomassRxnID))])
    disp(['                            after model reductionF' num2str(testSolx(1))])
end
if abs(solWTx(biomassRxnID)-testSolx(1)) >= 10^-5
    warning('Growth rate of wild type strain were incorrect.')
    reductionStatus = 2;
end

% FBA simulation on single knockout strains using original model;
rxnDelList = find(solWTx~=0);
tempFluxSol = cell(1,length(rxnDelList));
% if ~exist('matlabpool') || (matlabpool('size') == 0)
    for i = 1 : length(rxnDelList)
        modelDel = model;
        modelDel.lb(rxnDelList(i)) = 0;
        modelDel.ub(rxnDelList(i)) = 0;
        delSolx = solveModel(modelDel);
        tempFluxSol{i} = delSolx;
%         tempFluxSol{i} = glpk(-modelDel.c,modelDel.S,modelDel.b,modelDel.lb,modelDel.ub,repmat('S',length(modelDel.mets),1));
    end
% else
%     parfor i = 1 : length(rxnDelList)
%         modelDel = model;
%         modelDel.lb(rxnDelList(i)) = 0;
%         modelDel.ub(rxnDelList(i)) = 0;
%         delSol = optimizeCbModel(modelDel,'max');
%         tempFluxSol{i} = delSol.x;
% %         tempFluxSol{i} = glpk(-modelDel.c,modelDel.S,modelDel.b,modelDel.lb,modelDel.ub,repmat('S',length(modelDel.mets),1));
%     end
% end

emptyIDs = find(cellfun('isempty',tempFluxSol));
for i = 1: length(emptyIDs)
    tempFluxSol{emptyIDs(i)} = zeros(length(model.rxns),1);
end
fluxSol = cat(2,tempFluxSol{:});
fluxSol = round(fluxSol*10^6)/10^6;
fluxSol(isnan(fluxSol))=0;
fluxSol(:,fluxSol(biomassRxnID,:)<=model.lb(biomassRxnID)) = 0;
tempGrRate = fluxSol(biomassRxnID,:);
tempGrRatio = tempGrRate/solWTx(biomassRxnID);
grRatio = ones(length(model.rxns),1);
grRatio(solWTx~=0) = tempGrRatio;
clear fluxSolution tempFluxSolution tempGrRate tempGrRatio rxnDelList

% FBA simulation on single knockout strains using reduced model;
biomassRxnID = 1;
targetRxnID = 2;
oxygenRxnID = 3;
rxnDelListMove = find(testSolx~=0);
tempFluxSol = cell(1,length(rxnDelListMove));
% if ~exist('matlabpool') || (matlabpool('size') == 0)
    for i = 1 : length(rxnDelListMove)
        modelDel = modelMove;
        modelDel.lb(rxnDelListMove(i)) = 0;
        modelDel.ub(rxnDelListMove(i)) = 0;
        delSolx = solveModel(modelDel);
        tempFluxSol{i} = delSolx;
    end
% else
%     parfor i = 1 : length(rxnDelListMove)
%         modelDel = modelMove;
%         modelDel.lb(rxnDelListMove(i)) = 0;
%         modelDel.ub(rxnDelListMove(i)) = 0;
%         delSol = optimizeCbModel(modelDel,'max');
%         tempFluxSol{i} = delSol.x;
%     end
% end
emptyIDs = find(cellfun('isempty',tempFluxSol));
for i = 1: length(emptyIDs)
    tempFluxSol{emptyIDs(i)} = zeros(length(modelMove.rxns),1);
end
fluxSol = cat(2,tempFluxSol{:});
fluxSol(isnan(fluxSol))=0;
fluxSol = round(fluxSol*10^6)/10^6;
fluxSol(:,fluxSol(biomassRxnID,:)<=modelMove.lb(biomassRxnID)) = 0;
tempGrRate = fluxSol(biomassRxnID,:);
tempGrRatio = tempGrRate/testSolx(biomassRxnID);
grRatioMove = ones(length(modelMove.rxns),1);
grRatioMove(testSolx~=0) = tempGrRatio;
clear fluxSolution tempFluxSolution tempGrRate tempGrRatio rxnDelListMove

compGrRatio = ones(length(model.rxns),2);
compGrRatio(:,1) = grRatio;
for i = 1 : length(modelMove.rxns)
    compGrRatio(rxnAssocMat(i,:),2) = grRatioMove(i);
end
if sum((compGrRatio(:,1)-compGrRatio(:,2)).^2) < 10^5
    if options.verbFlag == true
        disp('Grwoth rates of all single knockout strains were correct.')
    end
else
    warning('Growth rates of some single knockout strain were incorrect.')
    reductionStatus = 3;
end
clear grRatio grRatioMove


%% Identify essential gene sets for either cell growth or target production

tempFluxSol1 = cell(1,length(geneSets));
tempFluxSol2 = cell(1,length(geneSets));
% if ~exist('matlabpool') || (matlabpool('size') == 0)
    for i = 1 : length(geneSets)
        modelDel = modelMove;
        modelDel.lb(logical(geneSetRxnMat(i,:))) = 0;
        modelDel.ub(logical(geneSetRxnMat(i,:))) = 0;
        delSolx = solveModel(modelDel);
        tempFluxSol1{i} = delSolx;
        modelDel = changeObjective(modelDel,modelDel.rxns(1:2),[0,1]);
        delSolx = solveModel(modelDel);
        tempFluxSol2{i} = delSolx;
    end
% else
%     parfor i = 1 : length(geneSets)
%         modelDel = modelMove;
%         modelDel.lb(logical(geneSetRxnMat(i,:))) = 0;
%         modelDel.ub(logical(geneSetRxnMat(i,:))) = 0;
%         delSol = optimizeCbModel(modelDel,'max');
%         tempFluxSol1{i} = delSol.x;
%         modelDel = changeObjective(modelDel,modelDel.rxns(1:2),[0,1]);
%         delSol = optimizeCbModel(modelDel,'max');
%         tempFluxSol2{i} = delSol.x;
%     end
% end
emptyIDs = find(cellfun('isempty',tempFluxSol1));
for i = 1: length(emptyIDs)
    tempFluxSol1{emptyIDs(i)} = zeros(length(modelMove.rxns),1);
end
emptyIDs = find(cellfun('isempty',tempFluxSol2));
for i = 1: length(emptyIDs)
    tempFluxSol2{emptyIDs(i)} = zeros(length(modelMove.rxns),1);
end
SKOfluxes1 = cat(2,tempFluxSol1{:});
SKOfluxes1(isnan(SKOfluxes1))=0;
SKOfluxes1 = round(SKOfluxes1*10^6)/10^6;
SKOfluxes2 = cat(2,tempFluxSol2{:});
SKOfluxes2(isnan(SKOfluxes2))=0;
SKOfluxes2 = round(SKOfluxes2*10^6)/10^6;
essentialGeneSets = (SKOfluxes1(1,:) <= modelMove.lb(biomassRxnID) | SKOfluxes2(2,:) == 0)';
clear tempFluxSolution delRxnSets


%% Arrange results of model reduction

modelReduced = modelMove;
modelReduced.rxnFormulas = rxnFormulas;
modelReduced.rxnAssociations = rxnAssociations;
modelReduced.rxnAssocMat = rxnAssocMat;
modelReduced.unqGeneSetsMat = unqGeneSetsMat;
modelReduced.geneSets = geneSets;
modelReduced.geneSetRxnMat = geneSetRxnMat;
modelReduced.geneSetAssocRxns = geneSetAssocRxns;
modelReduced.essentialGeneSets = essentialGeneSets;
modelReduced.oriRxns = model.rxns;
modelReduced.reductionStatus = reductionStatus;

biomassRxn = modelReduced.rxns{biomassRxnID};
targetRxn = modelReduced.rxns{targetRxnID};
oxygenRxn = modelReduced.rxns{oxygenRxnID};

if options.verbFlag == true
    if modelReduced.reductionStatus == 1
        disp('Model reduction was succeeded.')
    end
    toc
end

