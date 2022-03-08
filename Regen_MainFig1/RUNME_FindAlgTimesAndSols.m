% Note- sometimes, RobustKnock can cause errors when trying to initialise the
% COBRA toolbox as it can fail to initialise submodules correctly. To fix
% this, disable your computer's internet connection while running
% robustknock, as this will stop initCobraToolbox from trying to
% initialise/update submodules

clear all

% set 'toDo' to 'time' to run each algorithm until they discover a gc
% design, and then compare the time each algorithm took
% set 'toDo' to 'sols' to run each algorithm for 6h and compare the designs
% found

toDo='time'; % 'time'


saveData=true;

% parameters
nRepeats=3;
nDeletions=10;
modelName='iML1515.mat';
targetRxn='EX_succ_e';% 'EX_tyr__L_e' 'EX_pyr_e'
solver='gurobi';
forbiddenRxns={};%{'H2Ot','O2t','CO2t','EX_o2_e','O2tex','O2tpp','EX_co2_e','CO2tex','CO2tpp','EX_h2o_e','H2Otex','H2Otpp'};

% skip repeats of algorithms that are deterministic, as results shouldn't
% change between runs
skipDeterministics=false;

algorithmsToRun=[2,3,4,5,7];% {'OptKnock','RobustKnock','GCopt','FastPros','OptGene','Wayman','gcFront'}py
% OptKnock is now not used as strains are not guaranteed to be coupled (searches for maximum target synthesis at maximum growth, rather than minimum target synthesis), so isn't a valid comparison
% Wayman algorithm not used as it is not compatible with this model format, and also because it didn't work when I ran it- perhaps incompatible with this version of MATLAB?

timeLimit=6*60*60;


% either terminate when product synthesis exceeds 10^-6, or terminate when
% product synthesis is so high that it will never be reached
if isequal(toDo,'time')
    endTol=10^-6;
elseif isequal(toDo,'sols')
    endTol=1000000;
else
    error('Enter a valid identifier for toDo')
end



% load model
changeCobraSolver(solver,'all',0);
global CBTDIR
if contains(cd,'\')
    modelPath=[strrep(cd,'Regen_MainFig1','gcFront'),'\', modelName];
else
    modelPath=[strrep(cd,'Regen_MainFig1','gcFront'),'/', modelName];
end
model = readCbModel(modelPath);


% do not consider reactions that are not gene associated or are on the
% forbidden list
forbiddenInd=or(ismember(model.rxns,forbiddenRxns),ismember(model.grRules,''));


% getting experimentally validated essential genes
essGenes=fGetEssGenes;

% adding any gene that is predicted to be essential in silico to ess gene
% list

% find deletions that happen when gene is deleted
silicoEss=false(length(model.genes),1);
delMap=false(length(model.rxns),length(model.genes));
for a=1:length(model.genes)
    if ~ismember(model.genes(a),essGenes)
        rxnInds=find(model.rxnGeneMat(:,a)==1);
        if ~isempty(rxnInds)
            % gene catalyses at least one reaction in model
            x=true(size(silicoEss));
            x(a)=false;
            for b=1:length(rxnInds)
                delMap(rxnInds(b),a)=~eval(model.rules{rxnInds(b)});
            end
        end
    end
end

% test these deletions to find genes that are essential in silico
environment=getEnvironment();
parfor a=1:length(model.genes)
    % no deletions means either gene deletion won't disable any reactions, or
    % gene deletion was already on the essential gene list. Either way, no
    % need to test
    if sum(delMap(:,a))==0
        continue
    end
    
    restoreEnvironment(environment);
    m=model;
    m.lb(delMap(:,a))=0;
    m.ub(delMap(:,a))=0;
    s=optimizeCbModel(m);
    if s.f<10^-6
        silicoEss(a)=true;
    end
end
essGenes=[essGenes;model.genes(silicoEss)];

% removing reactions if they are the only reaction catalysed by an
% essential gene, as disruption of this reaction is likely what causes
% lethality after gene KO. If an essential gene affects multiple reactions,
% hard to know which of the reactions must be maintained
% for a=1:length(essGenes)
%     geneInd=find(ismember(model.genes,essGenes(a)));
%     if ~isempty(geneInd)
%         rxnInds=find(logical(model.rxnGeneMat(:,geneInd)));
%         if length(rxnInds)==1
%             forbiddenInd(rxnInds)=true;
%         end
%     end
% end

% removing reactions that cannot be KOd without deletion of an essential
% gene

for a=1:length(model.rxns)
    
    % remove reaction if flux is mandatory
    if model.lb(a)>0 || model.ub(a)<0
        forbiddenInd(a)=true;
        continue
    end
    
    % test if rxn is KO'd if every non-essential gene is KO'd- if it isn't,
    % then not a valid experimental target
    if forbiddenInd(a)==false
        x=ismember(model.genes,essGenes);
        cantKO=eval(model.rules{a});
        if cantKO
            % rxn can only be knocked out if an essential gene is knocked
            % out
            forbiddenInd(a)=true;
            
        else
            % test if rxn can be deleted
            m=model;
            m.lb(a)=0;
            m.ub(a)=0;
            s=optimizeCbModel(m);
            if s.f==0
                % rxn is necessary for growth- not a useful target
                forbiddenInd(a)=true;
            else
                % find non-essential genes affiliated with the reaction
                assocGenes=find(model.rxnGeneMat(a,:)==1);
                essentialAssocGenes=ismember(model.genes(assocGenes),essGenes);
                assocGenes=assocGenes(~essentialAssocGenes);
                
                % try and find a combination of KOs of non-essential genes
                % that KO the reaction without preventing growth
                forbiddenInd(a)=true;
                for b=length(assocGenes):-1:1
                    
                    % create all possible combinations of size b of the
                    % associated genes
                    combos=nchoosek(1:length(assocGenes),b);
                    for c=1:size(combos,1)
                        % knockout the genes
                        x=true(length(model.genes),1);
                        x(assocGenes(combos(c,:)))=false;
                        
                        % test effect of knocking out genes on reaction of
                        % interest
                        if eval(model.rules{a})
                            % combination doesn't disable reaction of
                            % interest
                            continue
                        else
                            % find other reactions affiliated with the
                            % genes
                            tempGeneInds=assocGenes(combos(c,:));
                            [tempRxnInds,~]=find(model.rxnGeneMat(:,tempGeneInds)==1);
                            
                            % see if other reactions are disabled if these
                            % genes are KOd
                            inactiveRxns=false(size(tempRxnInds));
                            for d=1:size(tempRxnInds)
                                inactiveRxns(d)=~eval(model.rules{tempRxnInds(d)});
                            end
                            tempRxnInds=tempRxnInds(inactiveRxns);
                            
                            % check if other reactions being disabled
                            % prevents growth
                            m=model;
                            m.lb(tempRxnInds)=0;
                            m.ub(tempRxnInds)=0;
                            s1=optimizeCbModel(m);
                            if s1.f>0
                                % a combination of gene KOs exists that will
                                % disable the reaction of interest without
                                % disabling growth
                                forbiddenInd(a)=false;
                                break
                            end
                        end
                    end
                    % stop loop once a combination that KOs the
                    % reaction without preventing growth is found
                    if forbiddenInd(a)==false
                        break
                    end
                end
                
            end
        end
    end
end




% deleting gene association for forbidden reactions to make sure FastPros won't use them
forbiddenInd=find(forbiddenInd);
for a=1:length(forbiddenInd)
    model.rules{forbiddenInd(a)}='';
    model.grRules{forbiddenInd(a)}='';
end
model.rxnGeneMat(forbiddenInd,:)=0;

forbiddenRxns=model.rxns(forbiddenInd);
allowedRxns=model.rxns(~ismember(model.rxns,forbiddenRxns));


origDirectory=cd;

nAlgorithms=7;

times=nan(nRepeats,nAlgorithms);
prodVals=zeros(nRepeats,nAlgorithms);

allDesigns=struct;


for a=1:nRepeats
    % close parallel pool- algorithms that use parallel processing should
    % include time needed to set this up
    delete(gcp('nocreate'));
    close all
    
    % shuffle order that algorithms are tested in- in case this makes a
    % difference
    algOrder=randperm(nAlgorithms);
    algOrder=algOrder(ismember(algOrder,algorithmsToRun));
    for b=1:length(algOrder)
        close all
        
        switch(algOrder(b))
            
            case 1
                % run tilted OptKnock- TILTING DOES NOT SEEM TO WORK
                % CORRECTLY
                fprintf('\nRunning OptKnock\n\n')
                cd([origDirectory,'/optknock'])
                startTime=tic;
                
                m=model;
                %                 m.c(ismember(m.rxns,targetRxn))=-0.001;
                %                 m.lb(ismember(m.rxns,targetRxn))=0.001;
                %
                
                optOK=struct;
                optOK.numDel=nDeletions;
                optOK.targetRxn=targetRxn;
                
                %                m.lb(m.c==1)=0.01;
                coptOK=struct;
                coptOK.rxnList={m.rxns(m.c==1)};
                coptOK.values=0.01;
                coptOK.sense='G';
                coptOK.rxnInd=find(m.c==1);
                
                solOK=OptKnock2(m,allowedRxns,optOK,coptOK,[],true);
                
                times(a,1)=toc(startTime);
                m2=model;
                if ~isempty(solOK.rxnList)
                    allDesigns(a).optknock=solOK.rxnList;
                    
                    m2.lb(ismember(m2.rxns,solOK.rxnList))=0;
                    m2.ub(ismember(m2.rxns,solOK.rxnList))=0;
                end
                m2.c(ismember(m2.rxns,targetRxn))=-0.001;
                s=optimizeCbModel(m2);
                prodVals(a,1)=s.x(ismember(m2.rxns,targetRxn));
                
                disp(solOK.rxnList')
                disp(prodVals(a,1))
                
                m2.c(ismember(m2.rxns,targetRxn))=0;
                prodEnv(m2,targetRxn);
                
                
            case 2
                % run RobustKnock- using OptPipe implementation as
                % RobustKnock implementation from original paper
                % requires TOMLAB
                if skipDeterministics && a>1
                    times(a,2)=times(1,2);
                    prodVals(a,2)=prodVals(1,2);
                else
                    
                    
                    fprintf('\nRunning RobustKnock\n\n')
                    cd([origDirectory,'/robustKnock/common_functions'])
                    
                    startTime=tic;
                    
                    [resultsRK,fluxSolutionsRK]=optPipe(model,model.rxns{model.c==1},targetRxn,nDeletions,'a',timeLimit,endTol, allowedRxns);
                    
                    times(a,2)=toc(startTime);
                    
                    if isempty(resultsRK)
                        disp('RobustKnock failed to find any solutions')
                    else
                        
                        if size(resultsRK,2)==1 && size(resultsRK,1)>1
                            % results output will be in wrong orientation
                            tRK=table(join(resultsRK',' '),fluxSolutionsRK(:,1),fluxSolutionsRK(:,2));
                        else
                            tRK=table(join(resultsRK,' '),fluxSolutionsRK(:,1),fluxSolutionsRK(:,2));
                        end
                        tRK=sortrows(tRK,3);
                        allDesigns(a).robustknock=tRK{end,1};
                        if tRK{end,3}>endTol
                            disp(tRK(end,:))
                            prodVals(a,2)=tRK{end,3};
                        else
                            disp('RobustKnock failed to find any solutions')
                        end
                        
                    end
                end
                
            case 3
                % run GCopt
                fprintf('\nRunning GCopt\n\n')
                cd([origDirectory,'/GCopt'])
                
                if skipDeterministics && a>1
                    times(a,3)=times(1,3);
                    prodVals(a,3)=prodVals(1,3);
                else
                    
                    
                    startTime=tic;
                    
                    modelGCO=model;
                    modelGCO.rev=modelGCO.lb<0;
                    %modelGCO.subSystems=[modelGCO.subSystems{:}]';
                    
                    probOpts.notKORxns  = model.rxns(~ismember(model.rxns,allowedRxns));
                    % assign biomass formation, substrate uptake and target reaction
                    probOpts.bmRxn              = modelGCO.rxns{modelGCO.c==1};
                    probOpts.subsRxn            = 'EX_glc__D_e';
                    probOpts.targetRxn          = targetRxn;
                    probOpts.sense  = 'max';
                    probOpts.maxKO      = nDeletions;
                    probOpts.fixMu      = [];
                    
                    genOpts.TimeLimit       = timeLimit;  % Maximal runtime of the solver [sec]
                    genOpts.compressFlag   = 1;
                    genOpts.modelType   = 0;
                    if strcmp(toDo,'time')
                        genOpts.BestObjStop= endTol; % added this term to allow termination after sufficient production achieved- is not exactly the same as the others, as GCopt assesses product synthesis level at the growth rate where max productivity occurs- but should be OK as long as long as the possible coupled designs have non-zero production at this value
                    else
                        genOpts.BestObjStop=inf;
                    end
                    
                    resultsGCO = init_gcOpt(modelGCO,probOpts,genOpts);
                    
                    times(a,3)=toc(startTime);
                    
                    if ~isempty(resultsGCO.KORxnNum)
                        modelGCO.lb(resultsGCO.KORxnNum)=0;
                        modelGCO.ub(resultsGCO.KORxnNum)=0;
                        modelGCO.c(ismember(modelGCO.rxns,targetRxn))=-0.001;
                        s=optimizeCbModel(modelGCO);
                        prodVals(a,3)=s.x(ismember(modelGCO.rxns,targetRxn));
                        disp(model.rxns(resultsGCO.KORxnNum)')
                        disp(prodVals(a,3))
                        
                        allDesigns(a).gcopt=model.rxns(resultsGCO.KORxnNum)';
                    end
                end
            case 4
                % run FastPros
                fprintf('\nRunning FastPros\n\n')
                cd([origDirectory,'/FastPros'])
                
                if skipDeterministics && a>1
                    times(a,4)=times(1,4);
                    prodVals(a,4)=prodVals(1,4);
                else
                    
                    m=model;
                    % editing grRules and rules so FP will consider rxn KOs and
                    % be comparable to other algorithms
                    geneAssoc=find(~ismember(m.grRules,''));
                    m.genes=cellstr(num2str([1:length(geneAssoc)]'));
                    m.rules=cellstr(strings(length(m.rxns),1));
                    m.grRules=cellstr(strings(length(m.rxns),1));
                    m.rxnGeneMat=m.rxnGeneMat*0;
                    for i=1:length(geneAssoc)
                        m.rules{geneAssoc(i)}=['x(',num2str(i),')'];
                        m.grRules{geneAssoc(i)}=num2str(i);
                        m.rxnGeneMat(geneAssoc(i),i)=1;
                    end
                    startTime=tic;
                    
                    m.rev=m.lb<0;
                    
                    
                    [modelFP,biomassRxnFP,targetRxnFP,oxygenRxnFP] = reduceModelForFP(m,m.rxns{m.c==1},targetRxn,'EX_o2_e');
                    
                    
                    
                    optFP=struct;
                    optFP.maxKONum=nDeletions;
                    optFP.timeLimit=timeLimit;
                    optFP.endTol=endTol;
                    optFP.rxnList=allowedRxns;
                    optFP.startTimeVar=startTime;
                    optFP.verbFlag=true;
                    
                    [solFP, statFP, resultFP] = FastPros2(modelFP, biomassRxnFP, targetRxnFP, oxygenRxnFP, optFP);
                    
                    times(a,4)=toc(startTime);
                    if isempty(solFP)
                        disp('FastPros did not find a solution')
                    else
                        disp(solFP(end).koRxnSets(1,:))
                        disp(solFP(end).prodRates(1))
                        prodVals(a,4)=solFP(end).prodRates(1);
                        
                        for c=1:size(solFP,2)
                            if c==1
                                allDesigns(a).fastpros=join(solFP(c).koRxnSets,' + ');
                            else
                                allDesigns(a).fastpros=[allDesigns(a).fastpros;join(solFP(c).koRxnSets,' + ')];
                            end
                        end
                    end
                end
            case 5
                % run OptGene (COBRA toolbox implementation)
                fprintf('\nRunning OptGene\n\n')
                cd([origDirectory,'/OptGene'])
                startTime=tic;
                
                [x, ~, ~, optGeneSol] = optGene2(model, targetRxn, 'EX_glc__D_e', allowedRxns, 'MaxKOs', nDeletions, 'TimeLimit', timeLimit, 'FitnessLimit', -endTol);
                
                times(a,5)=toc(startTime);
                if ~isfield(optGeneSol,'scores') || max(-optGeneSol.scores)<=0
                    disp('OptGene did not find a solution')
                else
                    disp(model.rxns(x)')
                    disp(max(-optGeneSol.scores))
                    prodVals(a,5)=max(-optGeneSol.scores);
                    
                    allDesigns(a).optgene=model.rxns(x)';
                end
            case 6
                % run Wayman- ALGORITHM DOESNT WORK! Ran into an error of
                % variable 'NA' not being set. Even if fixed, code seems to
                % require model in a completely different format, so can't make a
                % direct comparison.
                
                fprintf('\nRunning Wayman et al. algorithm\n\n')
                cd([origDirectory,'/Wayman'])
                startTime=tic;
                
                times(a,6)=toc(startTime);
            case 7
                % run gcFront (note that gcFront may be referred to as sgco
                % at points, as the algorithm name was changed after this script 
                % was written)
                
                fprintf('\nRunning gcFront\n\n')
                cd(strrep(origDirectory,'Regen_MainFig1','gcFront'))
                startTime=tic;
                
                optGCFRONT=struct;
                
                optGCFRONT.fitnesslimit=endTol;
                
                optGCFRONT.timelimit=timeLimit;
                optGCFRONT.maxknockouts=nDeletions;
                optGCFRONT.saveresults=false;
                optGCFRONT.newredundantremoval=false;
                
                optGCFRONT.ignorelistrxns=forbiddenRxns;
                optGCFRONT.dontkoess=false; % all computationally essential genes are in forbiddenRxns, and it would be an unfair comparison to add the time required to calculate these onto gcFront's time when this is not being added to the others
                [allTable,algParams,reducedModel]=gcFront(model,targetRxn,optGCFRONT);
                
                times(a,7)=toc(startTime);
                if isempty(allTable)
                    disp('gcFront did not find a solution')
                else
                    disp(max(allTable{:,4}))
                    prodVals(a,7)=max(allTable{:,4});
                    allDesigns(a).sgco=allTable{:,1};
                end
        end
        
        cd(origDirectory)
        
        save('tempdesigns.mat','allDesigns');
        save('tempprods.mat','prodVals');
        save('temptimes.mat','times');
        
        
    end
    
    
    
    if ~isfield(allDesigns(a),'robustknock')
        allDesigns(a).robustknock=[];
    end
    if ~isfield(allDesigns(a),'gcopt')
        allDesigns(a).gcopt=[];
    end
    if ~isfield(allDesigns(a),'optgene')
        allDesigns(a).optgene=[];
    end
    if ~isfield(allDesigns(a),'fastpros')
        allDesigns(a).fastpros=[];
    end
    if ~isfield(allDesigns(a),'sgco')
        allDesigns(a).sgco=[];
    end
       
    
end


if isequal(toDo,'sols')
    
    % find metrics for every identified design
    
    disp('Determining metrics for designs')
    
    cd(strrep(origDirectory,'Fig 1','gcFront'));
    
    targetInd=ismember(model.rxns,targetRxn);
    
    lpModel=convertModel(model);
    lpModel.c(targetInd)=-0.0001;
    
    poolModel=parallel.pool.Constant(lpModel);
    
    
    bigTable=table(" ",0,0,0,0,{" "},'VariableNames',{'Algorithm','Growth','Product','CouplingStrength','Deletions','Rxns'});
    bigTable(1,:)=[];
    
    for a=1:length(allDesigns)
        
        % finding metrics for robustknock
        if ~isempty(allDesigns(a).robustknock)
            x=ismember(model.rxns,splitString(allDesigns(a).robustknock{:},' '))';
            fitnesses=myFitness3obj(x, poolModel, 1:length(model.rxns), find(targetInd), find(model.c==1), length(model.rxns), 10^-8, 0, 10^-5, 0, 10^-4, false, [], []);
            fitnesses=-fitnesses;
            tempTable=table(repmat("RobustKnock",[size(x,1),1]),fitnesses(:,1),fitnesses(:,2),fitnesses(:,3),sum(x,2),join(model.rxns(x),' + '),'VariableNames',{'Algorithm','Growth','Product','CouplingStrength','Deletions','Rxns'});
            bigTable=[bigTable;tempTable];
        end
        
        % finding metrics for gcOpt
        if ~isempty(allDesigns(a).gcopt)
            x=ismember(model.rxns,allDesigns(a).gcopt)';
            fitnesses=myFitness3obj(x, poolModel, 1:length(model.rxns), find(targetInd), find(model.c==1), length(model.rxns), 10^-8, 0, 10^-5, 0, 10^-4, false, [], []);
            fitnesses=-fitnesses;
            tempTable=table(repmat("gcOpt",[size(x,1),1]),fitnesses(:,1),fitnesses(:,2),fitnesses(:,3),sum(x,2),join(model.rxns(x),' + '),'VariableNames',{'Algorithm','Growth','Product','CouplingStrength','Deletions','Rxns'});
            bigTable=[bigTable;tempTable];
        end
        
        % finding metrics for optGene
        if ~isempty(allDesigns(a).optgene)
            x=ismember(model.rxns,allDesigns(a).optgene)';
            fitnesses=myFitness3obj(x, poolModel, 1:length(model.rxns), find(targetInd), find(model.c==1), length(model.rxns), 10^-8, 0, 10^-5, 0, 10^-4, false, [], []);
            fitnesses=-fitnesses;
            tempTable=table(repmat("optGene",[size(x,1),1]),fitnesses(:,1),fitnesses(:,2),fitnesses(:,3),sum(x,2),join(model.rxns(x),' + '),'VariableNames',{'Algorithm','Growth','Product','CouplingStrength','Deletions','Rxns'});
            bigTable=[bigTable;tempTable];
        end
        
        % finding metrics for FastPros
        if ~isempty(allDesigns(a).fastpros)
            fpDels=strrep(allDesigns(a).fastpros, ',', '+');
            x=false(size(fpDels,1),length(model.rxns));
            for b=1:size(x,1)
                delString=splitString(fpDels{b},'+');
                delString=strrep(delString,' ','');
                for c=1:length(delString)
                    slashInd=regexp(delString{c},'/');
                    if ~isempty(slashInd)
                        delString{c}=delString{c}(1:slashInd-1);
                    end
                end
                x(b,:)=ismember(model.rxns,delString)';
            end
            
            fitnesses=myFitness3obj(x, poolModel, 1:length(model.rxns), find(targetInd), find(model.c==1), length(model.rxns), 10^-8, 0, 10^-5, 0, 10^-4, false, [], []);
            fitnesses=-fitnesses;
            tempTable=table(repmat("FastPros",[size(x,1),1]),fitnesses(:,1),fitnesses(:,2),fitnesses(:,3),sum(x,2),repmat({' '},[size(x,1),1]),'VariableNames',{'Algorithm','Growth','Product','CouplingStrength','Deletions','Rxns'});
            for b=1:size(tempTable,1)
                tempTable{b,6}=join(model.rxns(x(b,:)),' + ');
            end
            bigTable=[bigTable;tempTable];
        end
        
        % finding metrics for gcFront
        if ~isempty(allDesigns(a).sgco)
            sgcoDels=allDesigns(a).sgco;
            x=false(size(sgcoDels,1),length(model.rxns));
            for b=1:size(x,1)
                delString=splitString(char(sgcoDels(b)),' ');
                for c=1:length(delString)
                    slashInd=regexp(delString{c},'/');
                    if ~isempty(slashInd)
                        delString{c}=delString{c}(1:slashInd-1);
                    end
                end
                x(b,:)=ismember(model.rxns,delString)';
            end
            
            fitnesses=myFitness3obj(x, poolModel, 1:length(model.rxns), find(targetInd), find(model.c==1), length(model.rxns), 10^-8, 0, 10^-5, 0, 10^-4, false, [], []);
            fitnesses=-fitnesses;
            tempTable=table(repmat("gcFront",[size(x,1),1]),fitnesses(:,1),fitnesses(:,2),fitnesses(:,3),sum(x,2),repmat({' '},[size(x,1),1]),'VariableNames',{'Algorithm','Growth','Product','CouplingStrength','Deletions','Rxns'});
            for b=1:size(tempTable,1)
                tempTable{b,6}=join(model.rxns(x(b,:)),' + ');
            end
            bigTable=[bigTable;tempTable];
        end
        
    end
    
    
    wtEnv=prodEnvFast(model,targetRxn,0);
    
    cd(origDirectory)
    
    outputData=struct;
    outputData.bigTable=bigTable;
    outputData.times=times;
    outputData.targetRxn=targetRxn;
    outputData.nDeletions=nDeletions;
    outputData.forbiddenRxns=forbiddenRxns;
    outputData.timeLimit=timeLimit;
    outputData.algorithmsToRun=algorithmsToRun;
    outputData.model=model;
    
    if saveData==1
        filename=['designs ',targetRxn, ' ', strrep(modelName,'.mat',''),' ',char(datetime),'.mat'];
        filename=strrep(filename,':','_');
        save(filename,'outputData');
    end 
    
    bigTable=unique(bigTable);
    
    if sum(bigTable{:,3}==0)>0
        disp('Some designs have been identified by the algorithms but do not lead to coupling')
        disp(bigTable(bigTable{:,3}==0,:))
        bigTable=bigTable(bigTable{:,3}~=0,:);
    end
    
    
    bigTable=sortrows(bigTable,3);
    
    
    
    figure();
    
    axis square
    hold on
    plot(wtEnv(:,1),wtEnv(:,2),'k','LineWidth',3)
    scatter(bigTable{ismember(bigTable{:,1},'gcFront'),2},bigTable{ismember(bigTable{:,1},'gcFront'),3},200,[1,0,0],'Marker','.')
    scatter(bigTable{ismember(bigTable{:,1},'RobustKnock'),2},bigTable{ismember(bigTable{:,1},'RobustKnock'),3},200,[0,0.4,1],'Marker','square','LineWidth',2)
    scatter(bigTable{ismember(bigTable{:,1},'gcOpt'),2},bigTable{ismember(bigTable{:,1},'gcOpt'),3},200,[0,0.6,0.83],'Marker','diamond','LineWidth',2)
    scatter(bigTable{ismember(bigTable{:,1},'FastPros'),2},bigTable{ismember(bigTable{:,1},'FastPros'),3},200,[0,0.8,0.66],'Marker','pentagram','LineWidth',2)
    scatter(bigTable{ismember(bigTable{:,1},'OptGene'),2},bigTable{ismember(bigTable{:,1},'OptGene'),3},200,[0,1,0.5],'Marker','o','LineWidth',2)
    legend('WT','gcFront','RobustKnock','GCopt','FastPros','OptGene')
    xlabel('Growth (h^-^1)')
    ylabel('Product (mmol/g/h)')
    set(gcf,'Color','w')
    set(gca,'FontSize',16)
    hold off
    
    
else
    

    
    % % find stdev of algorithm run times
    stdevs=zeros(1,nAlgorithms);
    for a=1:nAlgorithms
        stdevs(a)=std(times(:,a)/60);
    end
    
    varNames={'OptKnock','RobustKnock','GCopt','FastPros','OptGene','Wayman','gcFront'};
    varNames=varNames(algorithmsToRun);
    
    % plot times with range as error bar
    f=figure();
    colormap cool
    hold on
    barNames=categorical(varNames);
    barNames=reordercats(barNames,varNames);
    for a=1:length(barNames)
        switch barNames(a)
            case 'OptKnock'
                bar(barNames(a), mean(times(:,1))/60);
            case 'RobustKnock'
                bar(barNames(a), mean(times(:,2)/60),'FaceColor', [0,0.4,1]);
            case 'GCopt'
                bar(barNames(a), mean(times(:,3)/60),'FaceColor', [0,0.6,0.83]);
            case 'FastPros'
                bar(barNames(a), mean(times(:,4)/60),'FaceColor', [0,0.8,0.66]);
            case 'OptGene'
                bar(barNames(a), mean(times(:,5)/60),'FaceColor', [0,1,0.5]);
            case 'Wayman'
                bar(barNames(a), mean(times(:,6))/60);
            case 'gcFront'
                bar(barNames(a), mean(times(:,7)/60),'FaceColor', [1,0,0]);
        end
    end
    
    er=errorbar(barNames,mean(times(:,algorithmsToRun)/60,1),stdevs(algorithmsToRun),'LineWidth',1);
    %er=errorbar(barNames,mean(times,1),mean(times,1)-min(times,[],1),max(times,[],1)-mean(times,1),'LineWidth',1);
    er.LineStyle='none';
    er.Color=[0,0,0];
    ylabel('Time (min)')
    set(f,'Color','w')
    set(gca,'FontSize',16)
    hold off
    
    prodTable=table(prodVals(:,1),prodVals(:,2),prodVals(:,3),prodVals(:,4),prodVals(:,5),prodVals(:,6),prodVals(:,7),'VariableNames',{'OptKnock','RobustKnock','GCopt','FastPros','OptGene','Wayman','gcFront'});
    fprintf('\n\n\n\nDisplaying minimum product synthesis for best design identified by each algorithm\n')
    disp(prodTable(:,algorithmsToRun))
    
    timeTable=table(times(:,1),times(:,2),times(:,3),times(:,4),times(:,5),times(:,6),times(:,7),'VariableNames',{'OptKnock','RobustKnock','GCopt','FastPros','OptGene','Wayman','gcFront'});
    disp('Displaying time taken by each algorithm')
    disp(timeTable(:,algorithmsToRun))
    
    outputData=struct;
    outputData.prodVals=prodVals;
    outputData.times=times;
    outputData.targetRxn=targetRxn;
    outputData.nDeletions=nDeletions;
    outputData.forbiddenRxns=forbiddenRxns;
    outputData.timeLimit=timeLimit;
    outputData.algorithmsToRun=algorithmsToRun;
    outputData.model=model;
    
    if saveData==1
        
        filename=['times ',targetRxn, ' ', strrep(modelName,'.mat',''),' ',char(datetime),'.mat'];
        filename=strrep(filename,':','_');
        save(filename,'outputData');
    end
    
end