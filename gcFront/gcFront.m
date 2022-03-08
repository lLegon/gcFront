function [designTable, options, model]=gcFront(model,targetRxn,options)
% [designTable, algParams, reducedModel]=gcFront(model,targetRxn,options)
% Uses multiobjective GA to find growth-coupled designs on the pareto front
% of growth, product synthesis and coupling strength
%
% INPUTS:
% model = either a structure containing a COBRA-toolbox compatible
% metabolic model, or the address of a file containing a COBRA-toolbox
% compatible metabolic model
%
% targetRxn = the reaction/metabolite in the model to be growth-coupled.
% if a metabolite name is supplied, then the exchange reaction for this
% metabolite will be selected
%
% options = a structure that contains parameters for the GA and
% pre/post-processing. Default values will be used for any fields that are
% not supplied by the user.
% options.solver = Character array of the name of the LP solver used for solving FBA problems. Default = currently set LP solver. If no solver set, will use whatever solver the function initCobraToolbox sets
% options.biomassrxn = Character array specifying the reaction that represents cell growth. Default = current objective reaction in model
% options.tol = Double specifying tolerance to mathematical errors. Flux differences smaller than this will be ignored. Default = 10^-8
% options.tiltval = Double specifying the coefficient used for objective value tilting. Default = 10^-4
% options.shiftval = Double specifying the size of the change in growth rate used to calculate shadow price of uncoupled designs. Default = 10^-5
% options.skipreduction = Logical that determines if algorithm should skip deletion of inactive reactions and pooling of unbranched pathways/redundant genes. Default = false (i.e. carry out reduction)
% options.mingrowth = Double specifying the minimum growth threshold- do not consider any deletion that lower growth below this value. Default = 10^-3
% options.minprod = Double specifying the minimum product threshold- do not consider any deletion that lowers maxmimum product synthesis below this value. Default = 10^-3
% options.removeredundancy = Logical that detemines if any redundant deletions should be removed from the designs that the GA identifies. Default = true
% options.saveresults = Logical that determines if results and algorithm parameters shoudl be saved. Default = true
% options.maxknockouts = Double specifying the maximum number of knockouts that a design may contain. Default = inf (i.e. no limit on the number of deletions allowed)
% options.deletegenes = Logical that determines if algorithm should test gene knockouts or reaction knockouts. Default = false (i.e. algorithm will knock out reactions)
% options.ignorelistrxns = cell array of reactions that user does not want to be knocked knocked out (e.g. reactions that are associated with experimentally inaccessible genes). Default = {} (i.e. no reactions)
% options.ignorelistgenes = cell array of genes that should not be knocked out (e.g. essential genes). If reaction knockouts are being tested, reactions that can't be knocked out without knocking out these genes will be ignored. Default = {} (i.e. no genes).
% options.dontkoess = Logical to determine whether reactions that can only be knocked out if an essential gene is knocked out are ignored. Default = true (i.e do not consider these reactions)
% options.onlykogeneassoc = Logical to determine if reactions that are not gene associated should be ignored. Default = true (i.e. only consider deletion of gene associated reactions)
% options.mutationrate = Double controlling the mean number of mutations that the GA mutation function introduces. Default = 1
% options.popsize = Double controlling the number of individuals in the GA population. Default = 200
% options.genlimit = Double specifying the maximum number of generations the GA will run for. Default = 10000
% options.timelimit = Double specifying the maximum length of time (in seconds) that the GA will run for. Default = 86400 (i.e. 1 day)
% options.fitnesslimit= Double. The algorithm will terminate if a design with a product synthesis rate higher than this is discovered. Default = inf (i.e. this condition will not stop the algorithm)
% options.spreadchangelimit = Double specifying how low the change in spread must be before the algorithm terminates. For further information see MATLAB's 'gamultiobj Algorithm' documentation. Default = 10^-4
% options.stallgenlimit = Double specifying the number of generations that change in spread must be below spreadchangelimit before the algorithm terminates. For further information see MATLAB's 'gamultiobj Algorithm' documentation. Default = number of generations (i.e. the algorithm will never terminate due to spread being too low)
% options.plotinterval = Double specifying the number of generations that must pass before the algorithm updates the plot of the current designs. Default = 1 (i.e. updates every generation)
% options.newredundantremoval = Logical specifying whether to use a new, faster method for removing redundant KOs, or to use the original version that was used to obtain the data for the gcFront paper.
% options.maxreductionsize = Double specifying the maximum size of design that should be fully explored for removal of redundant KOs
%
% OUTPUTS:
% designTable = a table containing deletions and metrics affiliated with
% those deletions
%
% algParams = the parameters that were used by the algorithm
%
% reducedModel = the model that was used by the GA. If the model was
% not reduced, this should be identical to the input model (except for
% changes to the objective)

% ensure that warnings are displayed
warning('on','all');

% checking global optimisation toolbox and COBRA toolbox are present
if isempty(which('gamultiobj.m'))
    error('MATLAB Global optimization toolbox is not available- please check that this package has been downloaded')
elseif isempty(which('initCobraToolbox.m'))
    error('COBRA toolbox is not available- please check that it has been downloaded and its path has been set')
end

% prompting user to enter a model and target if not supplied
if nargin==0 || and(isempty(model), isempty(targetRxn))
    promptAnswer=inputdlg({'Model file name/address (leave blank to select manually)','Target reaction/metabolite'},'gcFront',[1,50]);
    if isempty(promptAnswer)
        error('No target metabolite/reaction specified')
    end
    model=promptAnswer{1};
    targetRxn=promptAnswer{2};
elseif isempty(model)
    promptAnswer=inputdlg({'Model file name/address (leave blank to select manually)'},'gcFront',[1,50]);
    model=promptAnswer{1};
elseif nargin==1 || isempty(targetRxn)
    promptAnswer=inputdlg({'Target reaction/metabolite'},'gcFront',[1,50]);
    targetRxn=promptAnswer{1};
end

% prompting user to change options if necessary
if nargin<3
    
    % create an options structure and remove things that are not input options
    [model, targetRxn, options]=checkInputs(model, targetRxn, struct);
    options=rmfield(options,{'starttime','modelname','targetreaction'});
    
    % get the remaining fieldnames
    optionFields=sort(fieldnames(options));
    
    % change all fields into text
    for a=1:length(optionFields)
        
        % get the current value in the field
        currentVal=eval(['options.',optionFields{a},';']);
        
        % change the value into a char
        if ischar(currentVal)
            % is already text- no need to change
            continue
        elseif isnumeric(currentVal)
            eval(['options.',optionFields{a},' = num2str(options.',optionFields{a},');']);
        elseif islogical(currentVal)
            if currentVal
                eval(['options.',optionFields{a},' = num2str(1);']);
            else
                eval(['options.',optionFields{a},' = num2str(0);']);
            end
        elseif iscell(currentVal)
            % the only cell fields are ignorelistrxns and ignorelistgenes,
            % which will be empty at the moment as this is their default
            % value- so just set to an empty field
            eval(['options.',optionFields{a},' = char;']);
        end
        
    end
    
    % create the GUI for entering options
    set(0,'units','pixels');
    screensize= get(0,'screensize');
    optionFig=uifigure('Position',[mean(screensize([1,3]))-200,mean(screensize([2,4]))-110,400,220],'AutoResizeChildren','off','Name','gcFront options','CloseRequestFcn',@endOptionsEntry);
    optionFig.UserData=options;
    
    % populate GUI with objects
    
    uibutton(optionFig,'Position',[5,5,390,20],'Text','Start algorithm','ButtonPushedFcn',@(dd,~)endOptionsEntry(dd) );
    
    panel=uipanel(optionFig,'Position',[5,130,390,60]);
    lbl=uilabel(panel,'Position',[5 5 380 50],'Text',{'Name of the biomass reaction. Leave blank to use the current', 'model objective.'});
    
    textEntry=uitextarea(optionFig,'Position',[5,25,390,100]);
    textEntry.Value=options.biomassrxn;
    textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'biomassrxn');
    
    uidropdown(optionFig,'Items',optionFields,'Position',[5 195 390 20],'editable','on','ValueChangedFcn',@(dd,event)optionsDropdown(dd,lbl,optionFig,textEntry));
    
    % stop algorithm until user has finished entering options
    uiwait(optionFig)
    
    % store user inputs and close figure
    options=optionFig.UserData;
    delete(optionFig);
    
    % convert ignorelist reactions/genes into the appropriate format
    if isempty(options.ignorelistrxns)
        options=rmfield(options,'ignorelistrxns');
    else
        options.ignorelistrxns=char(strrep(options.ignorelistrxns,', ',','));
        options.ignorelistrxns=splitString(options.ignorelistrxns,',');
    end
    if isempty(options.ignorelistgenes)
        options=rmfield(options,'ignorelistgenes');
    else
        options.ignorelistgenes=char(strrep(options.ignorelistgenes,', ',','));
        options.ignorelistgenes=splitString(options.ignorelistgenes,',');
    end
    
end


% recording time that algorithm started
options.starttime=datetime('now','Format','yy_M_d HH_mm');

% checking that inputs are valid, and setting any unsupplied inputs to
% default values
[model, targetRxn, options]=checkInputs(model, targetRxn, options);

% dealing parameters in options to the appropriate variable

[tol, tiltVal, shiftVal, delOptions.skipReduction, minGrowth, minProd, removeRedundancy, saveResults, maxKnockouts, deleteGenes, delOptions.ignoreListRxns, delOptions.ignoreListGenes, delOptions.dontKoEss, delOptions.onlyKoGeneAssoc, popSize, genLimit, timeLimit, fitnessLimit, spreadChangeLimit, stallGenLimit, plotInterval, newredundantremoval, maxreductionsize]=deal( options.tol, options.tiltval, options.shiftval, options.skipreduction, options.mingrowth, options.minprod, options.removeredundancy, options.saveresults, options.maxknockouts, options.deletegenes, options.ignorelistrxns, options.ignorelistgenes, options.dontkoess, options.onlykogeneassoc, options.popsize, options.genlimit, options.timelimit, options.fitnesslimit, options.spreadchangelimit, options.stallgenlimit, options.plotinterval, options.newredundantremoval, options.maxreductionsize);

% reduce model

% store original model (for finding secretion later)
origModel=model;

[model,delsToInds]=quickReduceModel(model, model.rxns{model.c==1}, targetRxn, minGrowth, minProd, tol, deleteGenes, delOptions);

if isempty(delsToInds)
    error('There do not seem to be any possible deletions')
end

if ~delOptions.skipReduction
    disp('Model reduced')
    if deleteGenes
        disp(['Number of genes in reduced model: ',num2str(length(model.genes))])
        disp(['Number of genes removed: ',num2str(length(origModel.genes)-length(model.genes))])
        disp(['Number of genes that can be knocked out: ',num2str(length(delsToInds))])
    else
        disp(['Number of reactions in reduced model: ',num2str(length(model.rxns))])
        disp(['Number of reactions removed: ',num2str(length(origModel.rxns)-length(model.rxns))])
        disp(['Number of reactions that can be knocked out: ',num2str(length(delsToInds))])
    end
end


% scale mutation rate by number of knockouts so you average the desired
% number of mutations per iteration
mutationRate=options.mutationrate/length(delsToInds);


% create geneRules and geneMat variables
if deleteGenes
    geneRules=strrep(model.rules,' ','');
    geneMat=logical(full(model.rxnGeneMat));
    geneCounter=count(model.genes(delsToInds),'+')+1;
else
    geneRules=[];
    geneMat=[];
    geneCounter=[];
end

% getting index variables
targetInd=find(ismember(model.rxns,targetRxn));
growthInd=find(model.c==1);
nRxns=length(model.rxns);

% tilting model and converting into a parallel.pool.Constant object
lpModel=convertModel(model);
lpModel.c(targetInd)=-tiltVal;
poolModel=parallel.pool.Constant(lpModel);
% poolModel=lpModel

% finding production envelope
wtEnv=prodEnvFast(model,targetRxn,0);

% find fitness of the WT design
% defaultFitness=myFitness3obj(false(1,length(delsToInds)), poolModel, delsToInds, targetInd, growthInd, nRxns, tol, minGrowth, shiftVal, -inf, tiltVal, deleteGenes, geneRules, geneMat);
% defaultFitness=-defaultFitness(end);
defaultFitness=-1000000;


% create function handles for the GA
CreationFcn=@(delSetLength,FitnessFcn,options)myCreation(delSetLength,FitnessFcn,options,mutationRate, maxKnockouts, geneCounter);
FitnessFcn=@(x)myFitness3obj(x, poolModel, delsToInds, targetInd, growthInd, nRxns, tol, minGrowth, shiftVal, defaultFitness, tiltVal, deleteGenes, geneRules, geneMat);
CrossoverFcn=@(parents, options, genomeLength, ~, score, population)myCrossover(parents, options, genomeLength, [], score, population,maxKnockouts,geneCounter);
MutationFcn=@(parents,options,GenomeLength,FitnessFcn,state,score,population)myMutation(parents,options,GenomeLength,FitnessFcn,state,score,population,mutationRate,maxKnockouts,geneCounter);
PlotFcns={@(options,state,~,~)myPlotPareto3obj(defaultFitness,wtEnv,fitnessLimit,options,state,[],[])};

% create GA options structure
gaoptions=struct;
gaoptions.PopulationType='bitString';
gaoptions.PopulationSize=popSize;
gaoptions.Generations=genLimit;
gaoptions.TimeLimit=timeLimit;
gaoptions.StallGenLimit=stallGenLimit;
gaoptions.TolFun=spreadChangeLimit;
gaoptions.PlotInterval=plotInterval;
gaoptions.CreationFcn=CreationFcn;
gaoptions.CrossoverFcn=CrossoverFcn;
gaoptions.MutationFcn=MutationFcn;
gaoptions.PlotFcns=PlotFcns;
gaoptions.Vectorized='on';

% other GA options that we will not give users the option to set- default will always be used.
% if necessary they can be changed from here though

% gaoptions.StallTimeLimit
% gaoptions.InitialPopulation
% gaoptions.EliteCount
% gaoptions.ParetoFraction
% gaoptions.PopInitRange
% gaoptions.CrossoverFraction
% gaoptions.MigrationDirection
% gaoptions.MigrationInterval
% gaoptions.MigrationFraction
% gaoptions.StallTest
% gaoptions.TolCon
% gaoptions.InitialScores
% gaoptions.NonlinConAlgorithm
% gaoptions.InitialPenalty
% gaoptions.PenaltyFactor
% gaoptions.FitnessScalingFcn
% gaoptions.SelectionFcn
% gaoptions.DistanceMeasureFcn
% gaoptions.HybridFcn
% gaoptions.Display
% gaoptions.OutputFcns
% gaoptions.UseParallel


% start multiobjective GA
[x,~,~,~,population,score]=gamultiobj(FitnessFcn,length(delsToInds),[],[],[],[],[],[],[],gaoptions);


% get the metrics of the designs on the pareto
x=unique(x,'rows');

if ~deleteGenes
    designTable=table(strings(size(x,1),1),zeros(size(x,1),1),zeros(size(x,1),1),zeros(size(x,1),1),zeros(size(x,1),1),nan(size(x,1),1),'VariableNames',{'ReactionDeletions','NoOfDels','GrowthRate','ProductFlux','CouplingStrength','DistFromIdeal'});
else
    designTable=table(strings(size(x,1),1),zeros(size(x,1),1),zeros(size(x,1),1),zeros(size(x,1),1),zeros(size(x,1),1),nan(size(x,1),1),'VariableNames',{'GeneDeletions','NoOfDels','GrowthRate','ProductFlux','CouplingStrength','DistFromIdeal'});
end

% find metrics for any designs that were in final population
noMatchFound=false(size(designTable,1),1);
for a=1:size(designTable,1)
    if sum(x(a,:))==0
        designTable{a,1}="";
    elseif ~deleteGenes
        designTable{a,1}=join(model.rxns(delsToInds(x(a,:))),' ');
    else
        designTable{a,1}=join(model.genes(delsToInds(x(a,:))),' ');
    end
    if deleteGenes
        designTable{a,2}=sum(x(a,:))+count(designTable{a,1},'+');
    else
        designTable{a,2}=sum(x(a,:));
    end
    matchingPop=find(ismember(population,x(a,:),'rows'));
    if ~isempty(matchingPop)
        designTable{a,3:5}=-score(matchingPop(1),:);
    else
        noMatchFound(a)=true;
    end
end

% search for metrics just in case a pareto front design was not in the final population
% (in theory this shouldnt happen, but doesn't hurt to check)
noMatchFound=find(noMatchFound);
if ~isempty(noMatchFound)
    noMatchFitness=myFitness3obj(x(noMatchFound,:), poolModel, delsToInds, targetInd, growthInd, nRxns, tol, minGrowth, shiftVal, defaultFitness, tiltVal, deleteGenes, geneRules, geneMat);
    for a=1:length(noMatchFound)
        designTable{noMatchFound(a),3:5}=-noMatchFitness(a,:);
    end
end

if max(designTable{:,4})>tol
    if removeRedundancy==true
        disp('Removing redundant deletions from identified designs')
        % test removal of single genes/reactions from coupled designs to see if score is
        % maintained
        
        coupledInds=designTable{:,4}>tol;
        
        currDesigns=x(coupledInds,:);
        currFitnesses=designTable{coupledInds,3:5};
        
        if newredundantremoval
            designTable=newReduceDesigns(model, currDesigns, currFitnesses, FitnessFcn, tol, deleteGenes, delsToInds, maxreductionsize);
        else
            designTable=reduceDesigns(model, currDesigns, currFitnesses, FitnessFcn, tol, deleteGenes, delsToInds);
        end
        
        % remove designs with non-unique metrics if they have more deletions than
        % others
        
        redundantDesigns=false(size(designTable,1),1);
        for a=1:length(redundantDesigns)
            if redundantDesigns(a)==true
                continue
            end
            
            matchingMetrics=find(ismember(designTable{:,3:5} ,designTable{a,3:5}, 'rows'));
            
            if length(matchingMetrics)>1
                
                moreDels=designTable{matchingMetrics,2}~=min(designTable{matchingMetrics,2});
                redundantDesigns(matchingMetrics(moreDels))=true;
                
            end
            
        end
        
        designTable=designTable(~redundantDesigns,:);
    else
        coupledInds=designTable{:,4}>tol;
        designTable=designTable(coupledInds,:);
    end
end

% calculate Euclidian distance between normalised metrics of coupled designs and ideal point
% (where all metrics are at their maximum)
coupledDesigns=designTable{:,4}>tol;
normMetrics=designTable{coupledDesigns,3:5}./[max(wtEnv,[],1),2];
designTable{coupledDesigns,6}=sqrt( sum( (1-normMetrics).^2, 2) );


% replacing growthRate and ProductFlux with NaN if design is not coupled
if max(designTable{:,4}>tol)
    uncoupled=designTable{:,4}==0;
    designTable{uncoupled,[3,4]}=nan;
end


% sort table by product synthesis
designTable=sortrows(designTable,[-4,3,-5,2]);

% remove forbidden genes/reactions if they are part of a set, and separate
% out gene sets if every member must be deleted
if deleteGenes
    % separate out genes that are in sets where every member must be
    % deleted
    designTable{:,1}=strrep(designTable{:,1},'+',' ');
    
    % get all genes that are used in the obtained designs
    usedGeneNames=unique(splitString(char(strrep(join(designTable{:,1},' '),'/',' ')),' '));
    
    % compare these genes to the forbidden genes
    forbiddenUsedGenes=find(ismember(usedGeneNames,options.ignorelistgenes));
    for a=1:length(forbiddenUsedGenes)
        % the forbidden gene will only be occuring as part of a set of
        % options- so delete it from this set
        designTable{:,1}=strrep(designTable{:,1},['/',usedGeneNames{forbiddenUsedGenes(a)}],'');
        designTable{:,1}=strrep(designTable{:,1},[usedGeneNames{forbiddenUsedGenes(a)},'/'],'');
    end
else
    % get rxns used in the designs
    usedRxnNames=unique(splitString(char(strrep(join(designTable{:,1},' '),'/',' ')),' '));
    % compare these rxns to the forbidden rxns and to ones that arent gene
    % associated
    if options.onlykogeneassoc
        forbiddenUsedRxns=find( or( ismember(usedRxnNames,options.ignorelistrxns), ismember(usedRxnNames,origModel.rxns(ismember(origModel.rules,'')))  ) );
    else
        forbiddenUsedRxns=find(ismember(usedRxnNames,options.ignorelistrxns));
    end

    % Should probably also remove rxns from sets if they can only be KOd
    % with removal of computationally essential genes, or genes from
    % ignorelistgenes.
    
    % remove these reactions from the sets they appear in
    for a=1:length(forbiddenUsedRxns)
        % the forbidden gene will only be occuring as part of a set of
        % options- so delete it from this set
        designTable{:,1}=strrep(designTable{:,1},['/',usedRxnNames{forbiddenUsedRxns(a)}],'');
        designTable{:,1}=strrep(designTable{:,1},[usedRxnNames{forbiddenUsedRxns(a)},'/'],'');
    end
end


disp(designTable);


% update figure to only show pareto front
currPlot=findobj(get(gca,'Children'));
if max(designTable{:,4})>0
    % coupled designs are present
    
    colPoints=[0,0.01,1;
        0.55,0.8,0;
        0.3,0.9,0.6];
    
    colPoints2=[colPoints(1,1):(colPoints(1,2)-colPoints(1,1))/30:colPoints(1,2), colPoints(1,2):(colPoints(1,3)-colPoints(1,2))/30:colPoints(1,3);
        colPoints(2,1):(colPoints(2,2)-colPoints(2,1))/30:colPoints(2,2), colPoints(2,2):(colPoints(2,3)-colPoints(2,2))/30:colPoints(2,3);
        colPoints(3,1):(colPoints(3,2)-colPoints(3,1))/30:colPoints(3,2), colPoints(3,2):(colPoints(3,3)-colPoints(3,2))/30:colPoints(3,3);]';
    
    cols=colormap(colPoints2);
    
    colInds=ceil( designTable{:,5} * size(cols,1)/2 );
    uncoupledCols=colInds<=0;
    colInds(uncoupledCols)=1;
    
    currCols=cols(colInds,:);
    currCols(uncoupledCols,:)=1;
    
    currPlot(1).XData=designTable{:,3};
    currPlot(1).YData=designTable{:,4};
    currPlot(1).CData=currCols;
else
    % no coupled designs- just show best coupling strength that was
    % discovered
    currPlot(1).XData=designTable{:,3};
    currPlot(1).YData=designTable{:,5};
    currPlot(1).CData=zeros(size(designTable,1),3);
end

% enable dcm so designs can be clicked to show their metrics

dcm=datacursormode;
dcm.Enable='on';
dcm.DisplayStyle='window';

set(dcm,'UpdateFcn',{@showDesignInfo,designTable})

% remove existing buttons
currFig=gcf;
currFig.Name='Pareto Front';

buttonInds=false(size(currFig.Children));
for a=1:length(buttonInds)
    if isa(currFig.Children(a),'matlab.ui.control.UIControl')
        buttonInds(a)=true;
    end
end
delete(currFig.Children(buttonInds))

% enable dcm button

uicontrol('Parent',currFig,'Style','pushbutton','String','Toggle design clicking','Callback',@dcmFunction,'Position',[0 0 120 15]);


% add production envelope button

uicontrol('Parent',currFig,'Style','pushbutton','String','Prod Env','Callback',{@prodEnvButtonFunction,origModel,'',targetRxn,wtEnv,deleteGenes},'Position',[125 0 120 15]);


% find exchange reactions in original model
excRxns=findExcRxns(origModel);
% remove any that are not in reduced model (as these can't carry flux)
for a=1:length(excRxns)
    if excRxns(a)==true && sum(contains(model.rxns,origModel.rxns(a)))==0
        excRxns(a)=false;
    end
end
excRxns=find(excRxns);

% add show secretion button

uicontrol('Parent',currFig,'Style','pushbutton','String','Show secretion','Callback',{@showSecretionButtonFunction,origModel,'',deleteGenes,excRxns},'Position',[250 0 120 15]);


% changing name of any existing production envelope figures so they wont
% get caught up in any new production envelopes that are made
currFigures=findobj('type','figure');
for a=1:length(currFigures)
    if strcmp(currFigures(a).Name,'Production Envelopes')
        currFigures(a).Name=[currFigures(a).Name,' '];
    end
end

options.timerunfor=datetime('now','Format','yy_M_d HH_mm')-options.starttime;

% store the list of target kos
if deleteGenes
    options.validkos=model.genes(delsToInds);
else
    options.validkos=model.rxns(delsToInds);
end

if saveResults==1
    % save results and parameters in the current folder
    filename=[targetRxn,' ',char(options.starttime)];
    
    optionsText=evalc('disp(options)');
    fileID=fopen(['_Parameters ', filename,'.txt'],'w');
    fprintf(fileID,'%s',optionsText);
    fclose(fileID);
    
    writetable(designTable,['_Designs ',filename,'.csv']);
end

end