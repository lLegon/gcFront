function ...
    [tempResult, FastProsStatus, selKoStrainFlux, zeroFluxGeneSetIDs] ...
    = calcUtargetSelKoStrains(...
    model, modelGetUtarget, allowableKoGeneSetIDs, ...
    biomassRxnID, targetRxnID, options,...
    pre_selKoStrainFlux, pre_zeroFluxGeneSetIDs,...
    tempKoNum, FastProsResult, FastProsStatus)
%calcUtargetSelKoStrains  attempts to calculate utarget of knockout strains
%  and select knockout strains whose utarget is larget in strain sets.
%
% Modified by Laurence Legon (30/09/21) so code can be made to stop testing
% new knockouts after a coupled strain is found, or after a time limit has
% been exceeded.
% 
%     [tempResult, FastProsStatus, selKoStrainFlux, zeroFluxGeneSetIDs] ...
%     = calcUtargetSelKoStrains(...
%     model, modelGetUtarget, allowableKoGeneSetIDs, ...
%     biomassRxnID, targetRxnID, options,...
%     pre_selKoStrainFlux, pre_zeroFluxGeneSetIDs,...
%     tempKoNum, FastProsResult, FastProsStatus)
% 
% INPUTS
%  model                    COBRA model structure containing all required fields to perform FastPros.
%  modelGetUtarget          Model to calculate utarget
%  allowableKoGeneSetIDs    Gene set IDs of allowable knockouts
%  biomassRxnID             Reaction reresenting biomass objective function
%  targetRxnID              Reaction to be maximized by FastPros
%  options                  Optional inputs of FastPros 
%                           (Same structure to options in FastPros)
%  pre_selKoStrainFlux      Fluxes of previous knockout strains selected as parents in this generation
%  pre_zeroFluxGeneSetIDs   Gene set combinations whose knockout has no effect on previous results
%  tempKoNum                Current knockout number
%  FastProsResult           Current result of FastPros
%  FastProsStatus           Current status of FastPros
% 
% 
% OUTPUTS
%  tempResult               Temporary result of FastPros in this generation
%  FastProsStatus           Current status of FastPros
%  selKoStrainFlux          Fluxes of knockout strains selected as parents in the next generation
%  zeroFluxGeneSetIDs       Gene set combinations whose knockout has no effect on the current results
% 
% 
% Aug. 5th, 2013    Satoshi OHNO

%% Prepare variables

pre_tempResult = FastProsResult(tempKoNum);
pre_selKoGeneSetIDs = pre_tempResult.selKoGeneSetIDs;
pre_selKoStrainNum = size(pre_selKoGeneSetIDs,1);
pre_selKoStrainUtarget = pre_tempResult.selKoStrainUtarget;

varNames = {...
    'firstProdKoNum','firstUnqProdKoNum','maxProdStrain'...
    'existProdKoStrains','existEssentialKoStrains'};
for i = 1: length(varNames)
    eval([varNames{i} ' = FastProsStatus.' varNames{i} ';']);
end

selStrainNum = options.selStrainNum;
allowableKoGeneNum = length(allowableKoGeneSetIDs);


%% Identify gene set combinations whose knockout has no effect on previous results

if tempKoNum == 1
    hasNoEffects = ismember(allowableKoGeneSetIDs,pre_zeroFluxGeneSetIDs);
    koGeneSetIDs = allowableKoGeneSetIDs;
    growth = zeros(allowableKoGeneNum,1);
    growth(hasNoEffects) = -888;
    utarget = zeros(allowableKoGeneNum,1);
    utarget(hasNoEffects) = pre_selKoStrainUtarget;
    prodRates = zeros(allowableKoGeneNum,1);
else
    parKoGeneSetIDs = cell(pre_selKoStrainNum,1);
    parGrowth = cell(pre_selKoStrainNum,1);
    parUtarget = cell(pre_selKoStrainNum,1);
    parProdRates = cell(pre_selKoStrainNum,1);
    [uniqueKoGeneIDs,uniqueFirstIDs] = unique(pre_zeroFluxGeneSetIDs(:,1:tempKoNum-1),'rows','first');
    [~ ,uniqueLastIDs] = unique(pre_zeroFluxGeneSetIDs(:,1:tempKoNum-1),'rows','last');
    sliced_zeroFluxGeneSetIDs = pre_zeroFluxGeneSetIDs(:,tempKoNum);  
%     if ~exist('matlabpool') || (matlabpool('size') == 0)
        for i = 1 : pre_selKoStrainNum;
            tempKoGeneSets = zeros(allowableKoGeneNum,tempKoNum);
            tempKoGeneSets(:,1:tempKoNum-1) = repmat(pre_selKoGeneSetIDs(i,:),allowableKoGeneNum,1);
            tempKoGeneSets(:,tempKoNum) = allowableKoGeneSetIDs;
            tempGlowth = zeros(allowableKoGeneNum,1);
            tempUtarget = zeros(allowableKoGeneNum,tempKoNum);
            tempUtarget(:,1:tempKoNum-1) = repmat(pre_selKoStrainUtarget(i,:),allowableKoGeneNum,1);
            tempProdRates = zeros(allowableKoGeneNum,1);
            
            [~,intersectIDs]= intersect(uniqueKoGeneIDs,pre_selKoGeneSetIDs(i,:),'rows');
            hasNoEffects = ismember(allowableKoGeneSetIDs(:,1) , ...
                sliced_zeroFluxGeneSetIDs(uniqueFirstIDs(intersectIDs):uniqueLastIDs(intersectIDs)));
            tempGlowth(hasNoEffects) = -888;
            tempUtarget(hasNoEffects,tempKoNum) = tempUtarget(hasNoEffects,tempKoNum-1);
            
            parKoGeneSetIDs{i} = tempKoGeneSets;
            parGrowth{i} = tempGlowth;
            parUtarget{i} = tempUtarget;
            parProdRates{i} = tempProdRates;
        end
%     else
%         parfor i = 1 : pre_selKoStrainNum;
%             tempKoGeneSets = zeros(allowableKoGeneNum,tempKoNum);
%             tempKoGeneSets(:,1:tempKoNum-1) = repmat(pre_selKoGeneSetIDs(i,:),allowableKoGeneNum,1);
%             tempKoGeneSets(:,tempKoNum) = allowableKoGeneSetIDs;
%             tempGlowth = zeros(allowableKoGeneNum,1);
%             tempUtarget = zeros(allowableKoGeneNum,tempKoNum);
%             tempUtarget(:,1:tempKoNum-1) = repmat(pre_selKoStrainUtarget(i,:),allowableKoGeneNum,1);
%             tempProdRates = zeros(allowableKoGeneNum,1);
%             
%             [~,intersectIDs]= intersect(uniqueKoGeneIDs,pre_selKoGeneSetIDs(i,:),'rows');
%             hasNoEffects = ismember(allowableKoGeneSetIDs(:,1) , ...
%                 sliced_zeroFluxGeneSetIDs(uniqueFirstIDs(intersectIDs):uniqueLastIDs(intersectIDs)));
%             tempGlowth(hasNoEffects) = -888;
%             tempUtarget(hasNoEffects,tempKoNum) = tempUtarget(hasNoEffects,tempKoNum-1);
%             
%             parKoGeneSetIDs{i} = tempKoGeneSets;
%             parGrowth{i} = tempGlowth;
%             parUtarget{i} = tempUtarget;
%             parProdRates{i} = tempProdRates;
%         end
%     end
    koGeneSetIDs = cat(1,parKoGeneSetIDs{:});
    growth = cat(1,parGrowth{:});
    utarget = cat(1,parUtarget{:});
    prodRates = cat(1,parProdRates{:});
end
simulationList = find(growth == 0);


%% Identify gene set combinations whose knockout led target production in previous simulation

prodStrainKoNums = find(existProdKoStrains==1);
if ~isempty(prodStrainKoNums)
    for i = prodStrainKoNums
        prodKoGeneSetIDs = FastProsResult(i+1).productionKoStrains.koGeneSetIDs;
        for j = 1: size(prodKoGeneSetIDs,1)
            isProdStrains = ismember(koGeneSetIDs(simulationList,:),prodKoGeneSetIDs(j,1:i));
            growth(simulationList(sum(isProdStrains,2)>=i)) = -777;
            utarget(simulationList(sum(isProdStrains,2)>=i),end) = -777;
        end
    end
end
simulationList = find(growth == 0);


%% Identify essential gene set combinations for cell growth or target production

if tempKoNum >= 2
    essentialStrainKoNums = find(existEssentialKoStrains==1);
    if ~isempty(essentialStrainKoNums)
        for i = essentialStrainKoNums(1:end)
            pre_essentialKoStrains = FastProsResult(i+1).essentialGeneSetCombinations;
            par_essentialKoStrainIDs = cell(size(pre_essentialKoStrains,1),1);
            simulating_koGeneSetIDs = koGeneSetIDs(simulationList,:);
%             if ~exist('matlabpool') || (matlabpool('size') == 0)
                for j = 1: size(pre_essentialKoStrains,1)
                    hasNoEffects = ismember(simulating_koGeneSetIDs,pre_essentialKoStrains(j,:));
                    par_essentialKoStrainIDs{j} = find(sum(hasNoEffects,2)>=i);
                end
%             else
%                 parfor j = 1: size(pre_essentialKoStrains,1)
%                     hasNoEffects = ismember(simulating_koGeneSetIDs,pre_essentialKoStrains(j,:));
%                     par_essentialKoStrainIDs{j} = find(sum(hasNoEffects,2)>=i);
%                 end
%             end
            essentialKoStrainIDs = cat(1,par_essentialKoStrainIDs{:});
            growth(simulationList(essentialKoStrainIDs)) = -999;
            utarget(simulationList(essentialKoStrainIDs),end) = -999;
        end
    end
end
simulationList = find(growth == 0);


%% Find unique knockout combinations

if tempKoNum == 1
    preSimulationList = simulationList;
else
    sortKoGeneSetIDs = sort(koGeneSetIDs(simulationList,:),2);
    [uniqueKoGeneSetIDs,uniqueIDs]= unique(sortKoGeneSetIDs,'rows');
    uniqueMat = logical(sparse(length(simulationList),size(uniqueKoGeneSetIDs,1)));
    [~, uniqueLoc]= ismember(sortKoGeneSetIDs,uniqueKoGeneSetIDs,'rows');
    for i = 1 : length(uniqueLoc)
        uniqueMat(i,uniqueLoc(i)) = true; % column: original koGeneSetIDs, row: unique koGeneSetIDs
    end
    tempPreSimulationList = simulationList;
    simulationList = find(growth == 0);
    preSimulationList = simulationList;    
    simulationList = columnVector(intersect(tempPreSimulationList(uniqueIDs),preSimulationList));
    uniqueMat = uniqueMat(ismember(tempPreSimulationList,preSimulationList),:);
end

if isempty(simulationList)
    tempResult = struct([]);
    selKoStrainFlux = [];
    zeroFluxGeneSetIDs = [];
    if options.verbFlag == true
        disp('All knockout strains were already calcluated.')
    end
    return
end

%% Prepare variables for utarget caluculations

nonSimulationList = find(growth ~= 0);

flux = zeros(length(model.rxns),size(koGeneSetIDs,1));
if tempKoNum == 1
    flux(:,nonSimulationList) = repmat(pre_selKoStrainFlux,1,length(nonSimulationList));
else
    [~,sameFluxIDs] = ismember(koGeneSetIDs(nonSimulationList,1:tempKoNum-1) , pre_selKoGeneSetIDs,'rows');
    flux(:,nonSimulationList) = pre_selKoStrainFlux(:,sameFluxIDs);
end

maxKoNum = tempKoNum;


%% Caluculate utarget of single knockout strains

isKoGenes = sparse(length(simulationList),length(model.geneSets));
for i = 1 : length(simulationList)
    isKoGenes(i,koGeneSetIDs(simulationList(i),:)) = 1;
end
tempIsKoRxns = isKoGenes * model.geneSetRxnMat > 0;

parSol = cell(1,length(simulationList));

% if ~exist('matlabpool') || (matlabpool('size') == 0)
    for i = 1 : length(simulationList)
 
        % CODE ADDED TO ALLOW PREMATURE TERMINATION
        if toc(options.startTimeVar)>options.timeLimit
            % skip calculation and treat combination like it was infeasible
            % so the algorithm won't test anything new
            parSol{i}.x = zeros(length(modelGetUtarget.rxns),1);
            parSol{i}.f = 0;
            parSol{i}.y = zeros(length(modelGetUtarget.mets)+1,1);
            continue
        end

        modelDel = modelGetUtarget;
        modelDel.ub(tempIsKoRxns(i,:),1)=0;
        modelDel.lb(tempIsKoRxns(i,:),1)=0;
           
        [tempSolx, tempSolf, ~, tempSoly]= solveModel(modelDel);
        if tempSolf ~= 0
            parSol{i}.x = tempSolx;
            parSol{i}.f = tempSolf;
            parSol{i}.y = tempSoly;
        else
            parSol{i}.x = zeros(length(modelDel.rxns),1);
            parSol{i}.f = 0;
            parSol{i}.y = zeros(length(modelDel.mets)+1,1);
        end
    end
% else
%     parfor i = 1 : length(simulationList)
%         modelDel = modelGetUtarget;
%         modelDel.ub(tempIsKoRxns(i,:),1)=0;
%         modelDel.lb(tempIsKoRxns(i,:),1)=0;
%         tempSol= optimizeCbModel(modelDel,'max');
%         if tempSol.f ~= 0
%             parSol{i}.x = tempSol.x;
%             parSol{i}.f = tempSol.f;
%             parSol{i}.y = tempSol.y;
%         else
%             parSol{i}.x = zeros(length(modelDel.rxns),1);
%             parSol{i}.f = 0;
%             parSol{i}.y = zeros(length(modelDel.mets)+1,1);
%         end
%     end
% end
clear isKoGenes

% Convert unique results to the original nonunique resultsunique
if tempKoNum == 1
    sol = [parSol{:}];
else

        for i = 1 : length(simulationList)
            % CODE ADDED TO SPEED THIS UP IF TIME LIMIT EXCEEDED AS IT CAN TAKE DAYS TO DO THIS BIT
            % ONLY INCLUDE THE SOLUTION IF IT WILL FIT THE LATER CRITERIA
            % FOR A COUPLED STRAIN
            if toc(options.startTimeVar)<options.timeLimit || and(round(parSol{i}.y(end)*10^6)/10^6>=0,parSol{i}.x(biomassRxnID)>model.lb(biomassRxnID))
                % USE ORIGINAL CODE ONLY IF TIME IS NOT UP, OR STRAIN MEETS
                % THE CRITERIA FOR COUPLED STRAINS USED LATER IN THE CODE
                findID1 = uniqueMat(preSimulationList==simulationList(i),:);
                findID2 = uniqueMat(:,findID1);
                sol(findID2) = parSol{i};
            end
            
        end
        % CODE ADDED TO PUT AN INFEASIBLE RESULT INTO ANY SPACE IN SOL THAT
        % WAS SKIPPED OVER
        if toc(options.startTimeVar)>=options.timeLimit
            if ~exist('sol','var')
                % in case time ran out before sol was created
                sol=struct;
                sol.x=[];
                sol.f=[];
                sol.y=[];
            end
            for i=1:length(preSimulationList)
                % Treat combinations that were skipped over due to time
                % constraint as infeasible. Shouldn't change results, as
                % things that are skipped were either not tested or
                % uncoupled
                if size(sol,2)<i || isempty(sol(i).f)
                    sol(i).x= zeros(length(modelGetUtarget.rxns),1);
                    sol(i).f = 0;
                    sol(i).y = zeros(length(modelGetUtarget.mets)+1,1);
                end
            end
        end

end
clear parSol

preFlux = sparse(round(cat(2,[sol.x])*10^6)/10^6);
preFlux(isnan(preFlux))=0;
flux(:,preSimulationList) = preFlux;
objv = cat(2,[sol.f]);
objv(isnan(objv))=0;
shadowPrices = cat(2,[sol.y]);
shadowPrices(isnan(shadowPrices))=0;
shadowPrices = sparse(round(shadowPrices*10^6)/10^6);
objv(preFlux(biomassRxnID,:) <= model.lb(biomassRxnID)) = 0;
growth(preSimulationList)=objv';
utarget(preSimulationList,end) = shadowPrices(end,:)';
clear shadowPrices preFlux objv

% Create structure representing solution types of knockout strains
if tempKoNum ~= 1
    %Find knockout strains whose solutions are infeasible
    tempinfProdIDs = find(growth==0 | utarget(:,end)>=0);
    tempNoninfProdIDs = columnVector(setdiff(1:size(koGeneSetIDs,1),tempinfProdIDs)');
    
    %Find unique knockout strains
    [~,uniqueIDs] = unique(sort(koGeneSetIDs(tempinfProdIDs,:),2),'rows');
    koGeneSetIDs = [koGeneSetIDs(tempinfProdIDs(uniqueIDs),:);koGeneSetIDs(tempNoninfProdIDs,:)];
    growth = [growth(tempinfProdIDs(uniqueIDs));growth(tempNoninfProdIDs)];
    utarget = [utarget(tempinfProdIDs(uniqueIDs),:);utarget(tempNoninfProdIDs,:)];
    prodRates = [prodRates(tempinfProdIDs(uniqueIDs),:);prodRates(tempNoninfProdIDs,:)];
    flux = [flux(:,tempinfProdIDs(uniqueIDs)),flux(:,tempNoninfProdIDs)];
end

koStrainSolType.feasibleIDs = find(growth~=0);
koStrainSolType.infeasibleIDs = find(growth==0);
koStrainSolType.prodIDs = find(utarget(:,end)>=0 & growth>model.lb(biomassRxnID));
koStrainSolType.infProdIDs = columnVector(union(koStrainSolType.infeasibleIDs,koStrainSolType.prodIDs));
% CODE ADDED SO INFEASIBLE STRAINS ARE NOT ANALYSED- SHOULD SAVE TIME, AND
% SHOULDNT MATTER AS IT DOESNT MATTER WHICH COMBINATIONS ARE ESSENTIAL IF
% WE ARENT DOING ANOTHER ROUND ANYWAY
if toc(options.startTimeVar)>options.timeLimit
    koStrainSolType.infProdIDs=columnVector(koStrainSolType.prodIDs);
end
clear sol


%% Perform general FBA simulation on strains with infeasible solution or target production

if isempty(koStrainSolType.infProdIDs)
    tempProdStrainIDs = [];
    tempEssentialStrainIDs = [];
    essentialGeneSetCombinations = [];
else
    
    tempKoGeneSetIDs = koGeneSetIDs(koStrainSolType.infProdIDs,:);
    isKoGenes = zeros(length(koStrainSolType.infProdIDs),length(model.geneSets));
    for i = 1 : length(koStrainSolType.infProdIDs)
        isKoGenes(i,tempKoGeneSetIDs(i,:)) = 1;
    end
    tempIsKoRxns = isKoGenes * model.geneSetRxnMat > 0;
    parSolP = cell(1,length(koStrainSolType.infProdIDs));
%     if ~exist('matlabpool') || (matlabpool('size') == 0)

        for i = 1 : length(koStrainSolType.infProdIDs)
            
            modelProd = model;
            modelProd.ub(tempIsKoRxns(i,:),1)=0;
            modelProd.lb(tempIsKoRxns(i,:),1)=0;
            [tempSolx, tempSolf] = solveModel(modelProd);
            if tempSolf ~= 0
                parSolP{i}.x = tempSolx;
            else
                parSolP{i}.x = zeros(length(model.rxns),1);
            end
        end
%     else
%         parfor i = 1 : length(koStrainSolType.infProdIDs)
%             modelProd = model;
%             modelProd.ub(tempIsKoRxns(i,:),1)=0;
%             modelProd.lb(tempIsKoRxns(i,:),1)=0;
%             tempSol = optimizeCbModel(modelProd,'max');
%             if tempSol.f ~= 0
%                 parSolP{i}.x = tempSol.x;
%             else
%                 parSolP{i}.x = zeros(length(model.rxns),1);
%             end
%         end
%     end
    solP = [parSolP{:}];
    generalFbaFlux = cat(2,[solP.x]);
    generalFbaFlux(isnan(generalFbaFlux))=0;
    generalFbaFlux = round(generalFbaFlux*10^6)/10^6;
    generalFbaFlux = sparse(generalFbaFlux);
    
    %Idenntify strains whose knocked out gene set combination is essential to reach minimum cell growth
    koStrainSolType.essentialIDs = koStrainSolType.infProdIDs(generalFbaFlux(biomassRxnID,:) <= model.lb(biomassRxnID));

    %Idenntify strains whose knocked out gene set combination is essential for the target produciton
    koStrainSolType.imPossibleProdIDs = koStrainSolType.infProdIDs(generalFbaFlux(targetRxnID,:) == 0);
    
    %Integrate essential gene set combinations to reach minimun cell growth and to produce the target metabolite
    tempEssentialStrainIDs = columnVector(union(koStrainSolType.essentialIDs,koStrainSolType.imPossibleProdIDs));    
    if isempty(tempEssentialStrainIDs)
        essentialGeneSetCombinations = [];
    else
        essentialGeneSetCombinations = koGeneSetIDs(tempEssentialStrainIDs,:);
        existEssentialKoStrains(tempKoNum)=1;
        generalFbaFlux(:,generalFbaFlux(biomassRxnID,:) <= model.lb(biomassRxnID)) = 0;
    end
    
    %Identify strains in which some amount of target production is nessesary to reach the minimun cell growth
    isMustProdStrain = generalFbaFlux(biomassRxnID,:)'>model.lb(biomassRxnID) & growth(koStrainSolType.infProdIDs)==0 ...
        & generalFbaFlux(targetRxnID,:)'>0;
    koStrainSolType.mustProdIDs = koStrainSolType.infProdIDs(isMustProdStrain);
    if ~isempty(koStrainSolType.mustProdIDs)
        utarget(koStrainSolType.mustProdIDs,end) = inf;
        mustProdStrainFlux = generalFbaFlux(:,isMustProdStrain);
        prodRates(koStrainSolType.mustProdIDs) = generalFbaFlux(targetRxnID,isMustProdStrain)';
    else mustProdStrainFlux = [];
    end

    %Identify target produciton strains
    if ~isempty(koStrainSolType.prodIDs)
        [~,tempProIDs] = ismember(koStrainSolType.prodIDs,koStrainSolType.infProdIDs);
        prodStrainFlux = generalFbaFlux(:,tempProIDs);
        prodRates(koStrainSolType.prodIDs) = generalFbaFlux(targetRxnID,tempProIDs)';
    else prodStrainFlux = [];
    end
    
    tempProdStrainIDs = [koStrainSolType.prodIDs;koStrainSolType.mustProdIDs];
    prodStrainFlux = [prodStrainFlux,mustProdStrainFlux];
    
end
clear solP


%% Sort strains in the order in which produciton strains, non-production strains, and essential strains

tempUtargetStrainIDs = columnVector(setdiff((1:size(koGeneSetIDs,1))' , [tempProdStrainIDs;tempEssentialStrainIDs]));
tempProdStrainIDs = columnVector(tempProdStrainIDs);
tempUtargetStrainIDs = columnVector(tempUtargetStrainIDs);
tempEssentialStrainIDs = columnVector(tempEssentialStrainIDs);
sortIDs = [tempProdStrainIDs;...
    tempUtargetStrainIDs;...
    tempEssentialStrainIDs];

% Sort production strains in the order of the target production rates
[~,tempSortIDs] = sort(prodRates(tempProdStrainIDs) , 'descend');
sortIDs(1:length(tempProdStrainIDs)) = tempProdStrainIDs(tempSortIDs);
if ~isempty(koStrainSolType.infProdIDs)
    prodStrainFlux = prodStrainFlux(:,tempSortIDs);
end

% Sort non-production strains in the order of the utarget
if tempKoNum == 1 || options.selIncUtargetStrains == false
    [~,tempSortIDs] = sort(utarget(tempUtargetStrainIDs,end) , 'descend');
    sortIDs(length(tempProdStrainIDs)+(1:length(tempUtargetStrainIDs))) = ...
        tempUtargetStrainIDs(tempSortIDs);
else
    deltaUtarget = utarget(tempUtargetStrainIDs,end)-utarget(tempUtargetStrainIDs,end-1);
    increasedUtargetIDs = find(deltaUtarget>0 | utarget(tempUtargetStrainIDs,end)>=0);
    nonIncUtargetIDs = find(deltaUtarget<=0 & utarget(tempUtargetStrainIDs,end)<0);
    [~,tempSortIDs]= sort(utarget(tempUtargetStrainIDs(increasedUtargetIDs),end),'descend');
    sortIDs(length(tempProdStrainIDs)+(1:length(tempUtargetStrainIDs))) = ...
        [tempUtargetStrainIDs(increasedUtargetIDs(tempSortIDs));tempUtargetStrainIDs(nonIncUtargetIDs)];
end

% Sort strain information
koGeneSetIDs = koGeneSetIDs(sortIDs,:);
growth = growth(sortIDs);
utarget = utarget(sortIDs,:);
prodRates = prodRates(sortIDs,:);
flux = flux(:,sortIDs);
koStrainSolType.feasibleIDs = find(growth~=0);
koStrainSolType.infeasibleIDs = find(growth==0);
koStrainSolType.prodIDs = find(prodRates > 0);
if isempty(koStrainSolType.infProdIDs)
else
    koStrainSolType = rmfield(koStrainSolType,{'infProdIDs','essentialIDs','imPossibleProdIDs','mustProdIDs'});
end


% Remove knockout strain information of non-unique gene set combinations
if tempKoNum == 1
    selKoStrainIDs = columnVector(setdiff(koStrainSolType.feasibleIDs,koStrainSolType.prodIDs));
else    
    [~,uniqueIDs] = unique(sort(koGeneSetIDs,2),'rows','first');
    uniqueIDs = sort(uniqueIDs);
    koGeneSetIDs = koGeneSetIDs(uniqueIDs,:);
    growth = growth(uniqueIDs);
    utarget = utarget(uniqueIDs,:);
    prodRates = prodRates(uniqueIDs);
    flux = flux(:,uniqueIDs);
    koStrainSolType.feasibleIDs = find(growth~=0);
    koStrainSolType.infeasibleIDs = find(growth==0);
    koStrainSolType.prodIDs = find(prodRates > 0);
    
    tempSelKoStrainIDs = columnVector(setdiff(koStrainSolType.feasibleIDs,koStrainSolType.prodIDs));
    if options.selIncUtargetStrains == true
        selKoStrainIDs = ...
            tempSelKoStrainIDs(utarget(tempSelKoStrainIDs,end)-utarget(tempSelKoStrainIDs,end-1)>0);
    else selKoStrainIDs = tempSelKoStrainIDs(utarget(tempSelKoStrainIDs,end)>-500);
    end
end
if length(selKoStrainIDs) >= selStrainNum
    selKoStrainIDs = selKoStrainIDs(1:selStrainNum);
end

%% Extract information of production strains

prodStrainIDs = koStrainSolType.prodIDs;
if isempty(prodStrainIDs)
    productionKoStrains = struct([]);
else
    productionKoStrains.koNum = tempKoNum;
    productionKoStrains.koGeneSetIDs = koGeneSetIDs(prodStrainIDs,:);
    productionKoStrains.koGeneSets = model.geneSets(productionKoStrains.koGeneSetIDs);
    if size(productionKoStrains.koGeneSets,2) ~= tempKoNum
        productionKoStrains.koGeneSets = productionKoStrains.koGeneSets';
    end
    productionKoStrains.koRxnSets = findRxnSetsFromGeneSetIDs(model,productionKoStrains.koGeneSetIDs);
    productionKoStrains.utarget = utarget(prodStrainIDs,:);
    productionKoStrains.prodRates = prodRates(prodStrainIDs);
    productionKoStrains.flux = full(prodStrainFlux);
end


%% Identify gene set combinations in the next iteration whose knockout has no effect on the current result

if ~isempty(selKoStrainIDs)
    selKoStrainIDs = columnVector(selKoStrainIDs);
end
isKnockedOut=zeros(allowableKoGeneNum,size(selKoStrainIDs,1));
for i = 1 : allowableKoGeneNum
    tempRxnTFs = logical(model.geneSetRxnMat(allowableKoGeneSetIDs(i),:));
    isKnockedOut(i,:) = sum(flux(tempRxnTFs,selKoStrainIDs)~=0,1);
end
[noEffectsKoGeneSetIDs,sameStrainIDs] = find(isKnockedOut==0);
zeroFluxGeneSetIDs = [koGeneSetIDs(selKoStrainIDs(sameStrainIDs),:),allowableKoGeneSetIDs(noEffectsKoGeneSetIDs)]; 


%% Display production strains

if isempty(productionKoStrains)
    if options.verbFlag == true
        disp('No production strains were identified.')
    end
else
    maxID = 1;
    tempMaxProdRate = productionKoStrains.prodRates(maxID);
    if options.verbFlag == true
        disp([num2str(length(productionKoStrains.prodRates)) ' production strains were identified.'])
        disp(['Maximum target production rate is ' num2str(tempMaxProdRate) ' mmol/gDW/h '])
        disp(['utarget = ' num2str(productionKoStrains.utarget(maxID,tempKoNum))])
        for i = 1:tempKoNum
            disp([''  productionKoStrains.koGeneSets{maxID,i} ' gene set is knocked out.'])
            disp([' ' productionKoStrains.koRxnSets{maxID,i} ' reaction set is knocked out.'])
        end
    end
    existProdKoStrains(tempKoNum) = 1;
    if firstProdKoNum == 0
        firstProdKoNum = tempKoNum;
    end
    if firstUnqProdKoNum ==0
        if any(productionKoStrains.utarget(:,tempKoNum) > 0)
            firstUnqProdKoNum = tempKoNum;
        end
    end
    if tempMaxProdRate > maxProdStrain.prodRate
        maxProdStrain.koNum = tempKoNum;
        maxProdStrain.koGeneSetIDs = productionKoStrains.koGeneSetIDs(1,:);
        maxProdStrain.koGeneSets = productionKoStrains.koGeneSets(1,:);
        maxProdStrain.koRxnSets = productionKoStrains.koRxnSets(1,:);
        maxProdStrain.prodRate = tempMaxProdRate;
        maxProdStrain.utarget = productionKoStrains.utarget(1,:);
        maxProdStrain.flux = productionKoStrains.flux(:,1);
    end
end

if options.verbFlag == true
    if tempKoNum == 1 || length(selKoStrainIDs) < selStrainNum
        disp(['All strains (' num2str(length(selKoStrainIDs)) ') go to the next generation.'])
    else
        if options.selIncUtargetStrains == true
            selNum = length(find(utarget(:,end)<0 & utarget(:,end)-utarget(:,end-1) > 0));
            disp([num2str(selStrainNum) '/' num2str(selNum) ...
                ' (' num2str(selStrainNum/selNum*100,2) '%) strains go to the next generation.' ])
        else
            selNum = length(find(utarget(:,end)<0> 0));
            disp([num2str(selStrainNum) '/' num2str(selNum) ...
                ' (' num2str(selStrainNum/selNum*100,2) '%) strains go to the next generation.' ])
        end
    end
end


%% Arrange variables

selKoGeneSetIDs = koGeneSetIDs(selKoStrainIDs,:);
selKoStrainUtarget = utarget(selKoStrainIDs,:);
selKoStrainFlux = flux(:,selKoStrainIDs);

tempResult.koNum = tempKoNum;
tempResult.selKoGeneSetIDs = selKoGeneSetIDs;
tempResult.selKoStrainUtarget = selKoStrainUtarget;
tempResult.productionKoStrains = productionKoStrains;
tempResult.essentialGeneSetCombinations = essentialGeneSetCombinations;

FastProsStatus.maxKoNum = maxKoNum;
FastProsStatus.firstProdKoNum = firstProdKoNum;
FastProsStatus.firstUnqProdKoNum = firstUnqProdKoNum;
FastProsStatus.maxProdStrain = maxProdStrain;
FastProsStatus.existProdKoStrains = existProdKoStrains;
FastProsStatus.existEssentialKoStrains = existEssentialKoStrains;
