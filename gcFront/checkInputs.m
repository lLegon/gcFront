function [model, targetRxn, options]=checkInputs(model, targetRxn, options)
% [model, targetRxn, options]=checkInputs(model, targetRxn, options)
% Checks the inputs to gcFront. If inputs are invalid, then
% warnings will be supplied and default value used
% 
% model- either a COBRA toolbox compatible model, or a file name/address for a COBRA toolbox compatible model
% targetRxn- the reaction to be growth coupled. If a metabolite name is supplied, then an exchange reaction involving this metabolite will be set as the target
% options- a structure where the fields are algorithm parameters, and the values are the desired parameter values

%% CHECKING OPTIONS

if ~isstruct(options)
    error('Options should be supplied in the form of a struct')
end

optionsnames=fieldnames(options);

validOptions={'solver' 'biomassrxn' 'tol' 'tiltval' 'shiftval' 'skipreduction' 'mingrowth' 'minprod' 'removeredundancy' 'saveresults' 'maxknockouts' 'deletegenes' 'ignorelistrxns' 'ignorelistgenes' 'dontkoess' 'onlykogeneassoc' 'mutationrate' 'popsize' 'genlimit' 'timelimit' 'fitnesslimit' 'spreadchangelimit' 'stallgenlimit' 'plotinterval' 'starttime','newredundantremoval','maxreductionsize'};

% check for any fields in the options structure that are not recognised
unknownOptions=optionsnames(~ismember(optionsnames,validOptions));
for a=1:length(unknownOptions)
    if ismember(lower(unknownOptions{a}),validOptions)
        if ~isfield(options,lower(unknownOptions{a}))
            eval(['options.',lower(unknownOptions{a}),'=options.',unknownOptions{a},';']);
        else
            warning(['Overlapping fields ',unknownOptions{a},' and ', lower(unknownOptions{a}),' exist in the options structure- using value from ',lower(unknownOptions{a})])
        end
    else
        warning(['Unrecognised option: ', unknownOptions{a},' was used'])
    end
end

% reset optionsnames with any new fields that have been added
optionsnames=fieldnames(options);

% create a blank field for any option that was not supplied
suppliedFields=validOptions(ismember(validOptions,optionsnames));
missingFields=validOptions(~ismember(validOptions,optionsnames));
for a=1:length(missingFields)
    eval(['options.',missingFields{a},' = [];']);
end


% check each supplied field is in the correct format, and put empty vectors in place of invalid values (which will later be replaced by default value)
for a=1:length(suppliedFields)
    switch suppliedFields{a}
        case 'solver'
            options.solver=checkInputType(options.solver,'char',suppliedFields{a},[1,inf]);
        case 'biomassrxn'
            options.biomassrxn=checkInputType(options.biomassrxn,'char',suppliedFields{a},[1,inf]);
        case 'tol'
            options.tol=checkInputType(options.tol,'posNumeric',suppliedFields{a},[1,1]);
        case 'tiltval'
            options.tiltval=checkInputType(options.tiltval,'posNumeric',suppliedFields{a},[1,1]);
        case 'shiftval'
            options.shiftval=checkInputType(options.shiftval,'posNumeric',suppliedFields{a},[1,1]);
        case 'skipreduction'
            options.skipreduction=checkInputType(options.skipreduction,'logical',suppliedFields{a},[1,1]);
        case 'mingrowth'
            options.mingrowth=checkInputType(options.mingrowth,'posNumeric',suppliedFields{a},[1,1]);
        case 'minprod'
            options.minprod=checkInputType(options.minprod,'numeric',suppliedFields{a},[1,1]);
        case 'removeredundancy'
            options.removeredundancy=checkInputType(options.removeredundancy,'logical',suppliedFields{a},[1,1]);
        case 'saveresults'
            options.saveresults=checkInputType(options.saveresults,'logical',suppliedFields{a},[1,1]);
        case 'maxknockouts'
            options.maxknockouts=checkInputType(options.maxknockouts,'posInteger',suppliedFields{a},[1,1]);
        case 'deletegenes'
            options.deletegenes=checkInputType(options.deletegenes,'logical',suppliedFields{a},[1,1]);
        case 'ignorelistrxns'
            options.ignorelistrxns=checkInputType(options.ignorelistrxns,'cell',suppliedFields{a},[1,inf]);
        case 'ignorelistgenes'
            options.ignorelistgenes=checkInputType(options.ignorelistgenes,'cell',suppliedFields{a},[1,inf]);
        case 'dontkoess'
            options.dontkoess=checkInputType(options.dontkoess,'logical',suppliedFields{a},[1,1]);
        case 'onlykogeneassoc'
            options.onlykogeneassoc=checkInputType(options.onlykogeneassoc,'logical',suppliedFields{a},[1,1]);
        case 'mutationrate'
            options.mutationrate=checkInputType(options.mutationrate,'posNumeric',suppliedFields{a},[1,1]);
        case 'popsize'
            options.popsize=checkInputType(options.popsize,'posInteger',suppliedFields{a},[1,1]);
        case 'genlimit'
            options.genlimit=checkInputType(options.genlimit,'posInteger',suppliedFields{a},[1,1]);
        case 'timelimit'
            options.timelimit=checkInputType(options.timelimit,'posNumeric',suppliedFields{a},[1,1]);
        case 'fitnesslimit'
            options.fitnesslimit=checkInputType(options.fitnesslimit,'numeric',suppliedFields{a},[1,1]);
        case 'spreadchangelimit'
            options.spreadchangelimit=checkInputType(options.spreadchangelimit,'posNumeric',suppliedFields{a},[1,1]);
        case 'stallgenlimit'
            options.stallgenlimit=checkInputType(options.stallgenlimit,'posInteger',suppliedFields{a},[1,1]);
        case 'plotinterval'
            options.plotinterval=checkInputType(options.plotinterval,'posInteger',suppliedFields{a},[1,1]);
        case 'newredundantremoval'
            options.newredundantremoval=checkInputType(options.newredundantremoval,'logical',suppliedFields{a},[1,1]);
        case 'maxreductionsize'
            options.maxreductionsize=checkInputType(options.maxreductionsize,'posNumeric',suppliedFields{a},[1,1]);
    end
end

% set empty options fields to defaults

if isempty(options.solver)
    % find if a solver has already been set for LP- if so, continue using this solver
    solverInfo=parseSolverParameters('LP');
    
    % if no solver is loaded, initialise cobra toolbox and use default solver
    if isempty(solverInfo.solver)
        disp('No solver detected- will initialise COBRA toolbox and then use the default LP solver')
        initCobraToolbox(false);
        solverInfo=parseSolverParameters('LP');
        if isempty(solverInfo.solver)
            error('COBRA toolbox failed to set a LP solver- please supply a valid solver name')
        end
    end
    
    % store the solver in the options structure
    options.solver=solverInfo.solver;
else
    
    % set to supplied solver
    solverResult=changeCobraSolver(options.solver,'LP',0);
    
    % if solver didn't work, show the error message
    if solverResult~=1
        changeCobraSolver(options.solver,'LP');
        error('Something went wrong while setting solver')
    end
    
end
if isequal(options.solver,'glpk')
    warning('Solver glpk may not work correctly with parallel processing- please consider using a different solver')
end
% options.biomassrxn will be dealt with after model is loaded
if isempty(options.tol)
    options.tol=10^-8;
end
if isempty(options.tiltval)
    options.tiltval=10^-4;
end
if isempty(options.shiftval)
    options.shiftval=10^-5;
end
if isempty(options.skipreduction)
    options.skipreduction=false;
end
if isempty(options.mingrowth)
    options.mingrowth=10^-3;
end
if isempty(options.minprod)
    options.minprod=10^-2;
end
if isempty(options.removeredundancy)
    options.removeredundancy=true;
end
if isempty(options.saveresults)
    options.saveresults=true;
end
if isempty(options.maxknockouts)
    options.maxknockouts=inf;
end
if isempty(options.deletegenes)
    options.deletegenes=false;
end
if isempty(options.ignorelistrxns)
    options.ignorelistrxns={};
end
if isempty(options.ignorelistgenes)
    options.ignorelistgenes={};
end
if isempty(options.dontkoess)
    options.dontkoess=true;
end
if isempty(options.onlykogeneassoc)
    options.onlykogeneassoc=true;
end
if isempty(options.mutationrate)
    options.mutationrate=1;
end
if isempty(options.popsize)
    options.popsize=200;
end
if isempty(options.genlimit)
    options.genlimit=10000;
end
if isempty(options.timelimit)
    options.timelimit=60*60*24;
end
if isempty(options.fitnesslimit)
    options.fitnesslimit=inf;
end
if isempty(options.spreadchangelimit)
    options.spreadchangelimit=10^-4;
end
if isempty(options.stallgenlimit)
    options.stallgenlimit=options.genlimit;
end
if isempty(options.plotinterval)
    options.plotinterval=1;
end
if isempty(options.newredundantremoval)
    options.newredundantremoval=true;
end
if isempty(options.maxreductionsize)
    options.maxreductionsize=15;
elseif options.newredundantremoval==false
    if options.maxreductionsize~=15 && options.removeredundancy==true
        warning('Old removal of redundant deletions has been selected, but this cannot use maxreductionsize.')
    end
    options.maxreductionsize=inf;
end
if options.maxreductionsize<options.maxknockouts && options.removeredundancy==true && ~isempty(optionsnames)
    warning(['Maximum number of allowable knockouts (',num2str(options.maxknockouts) ,') is greater than maximum number of knockouts than are allowed during design reduction (',num2str(options.maxreductionsize),'), so some designs may not have redundant KOs removed'])
end


%% CHECKING MODEL

global CBTDIR

% if supplied model is a link to a file, then load model from this file
if iscell(model)
    model=model{:};
end
if isstring(model)
    model=char(model);
end
if ischar(model)
    try
        % first try loading model from current directory
        if ~isempty(model)
            if contains(cd,'/')
                model=readCbModel([cd,'/',model]);
            else
                model=readCbModel([cd,'\',model]);
            end
        else
            error('Input must be a row vector of characters')
        end
    catch errorStruct
        % if nothing found, just try loading model without specifying where
        try
            model=readCbModel(model);
        catch errorStruct2
            if length(model)>=4 && isequal(model(end-3:end),'.mat')
                % if still nothing and model is a .mat file, try extracting
                % only the necessary fields from this file in case
                % some redundant field is causing problems
                try
                    model=load(model);
                    tempFields=fields(model);
                    if length(tempFields)~=1
                        error(' ')
                    else
                        eval(['model=model.',tempFields{:},';'])
                    end
                    necessaryFields={'S','mets','b','csense','rxns','lb','ub','c','osenseStr','genes','rules','grRules','rxnGeneMat','modelID'};
                    redundantFields=setdiff(fields(model),necessaryFields);
                    if ~isempty(redundantFields)
                        model=rmfield(model,redundantFields);
                    end
                catch
                    % terminate if model could not be loaded
                    warning(['Searched current directory for model- error message: ',errorStruct.message])
                    warning(['Searched provided address for model- error message: ',errorStruct2.message])
                    error('Model address does not lead to a valid model')
                end
                
            else
                % terminate if model could not be loaded
                warning(['Searched current directory for model- error message: ',errorStruct.message])
                warning(['Searched provided address for model- error message: ',errorStruct2.message])
                error('Model address does not lead to a valid model')
            end
            
        end
    end
end

if isstruct(model)
    % check that model is valid, and if it is not, then check that the
    % reason for failure wasn't that the model is missing a genes field
    verifyResults=verifyModel(model);
    if ~isempty(fieldnames(verifyResults)) && sum(ismember(verifyResults.Errors.missingFields,{'genes','rules'}))~=length(verifyResults.Errors.missingFields)
        error('supplied model is a struct, but not a valid model')
    end
else
    % if model is not a structure, then it is not valid
    error('Model should either be a struct, or a character array specifying the filepath to a model that can be read with readCbModel')
end

% if biomass reaction was supplied, check that it is in model
if ~isempty(options.biomassrxn) && ~ismember(options.biomassrxn,model.rxns)
    warning(['Supplied biomass reaction: ',options.biomassrxn,' is not in the model- default biomass reaction will be used'])
    options.biomassrxn=[];
end

if isempty(options.biomassrxn)
    % if biomass reaction not supplied, return an error if no default objective found
    biomassInd=model.c==max(model.c);
    if max(model.c)<=0
        error('No default biomass reaction detected- please specify a biomass reaction')
    elseif sum(biomassInd)>1
        error('Multiple biomass reactions detected- please specify a biomass reaction')
    end
    
    % if no biomass reaction specified, then use the default objective as the biomass reaction
    options.biomassrxn=model.rxns{biomassInd};
end
model=changeObjective(model,options.biomassrxn);

% check if the fields rxnGeneMat, genes, rules and grRules are present
geneFields={'rxnGeneMat','genes','rules','grRules'};
presentGeneFields=isfield(model,geneFields);
if sum(presentGeneFields)~=length(presentGeneFields)
    % if they are not, check if any options that use these fields are
    % enabled
    if options.onlykogeneassoc || options.deletegenes || ~isempty(options.ignorelistgenes) || options.dontkoess
        
        if isfield(model,'genes') && or(isfield(model,'rules'),isfield(model,'grRules'))
            % if model contains genes and at least 1 rules field, build the other
            % fields from this

            if ~isfield(model,'grRules')
                model.grRules=model.rules;
                
                model.grRules=strrep(model.grRules,'&','and');
                model.grRules=strrep(model.grRules,'|','or');
                
                for a=1:length(model.genes)
                    model.grRules=strrep(model.grRules,['x(',num2str(a),')'],model.genes{a});
                end
            end
            
            if ~isfield(model,'rules')
                model.rules=model.grRules;
                model.rules=strrep(model.rules,'and','&');
                model.rules=strrep(model.rules,'or','|');
                
                for a=1:length(model.genes)
                    model.rules=strrep(model.rules,model.genes{a},['x(',num2str(a),')']);
                end
            end
            
            if ~isfield(model,'rxnGeneMat')
                model=buildRxnGeneMat(model);
            end
                
        else
            
            % warn the user if this will change a setting that they supplied
            geneUsingFields={'onlykogeneassoc','deletegenes','ignorelistgenes','dontkoess'};
            userSetFields=ismember(geneUsingFields,suppliedFields);
            if sum(userSetFields)~=0
                warning(['Model does not contain field(s): ',geneFields{~presentGeneFields}, ' so option(s): ',geneUsingFields{userSetFields},' cannot be used'])
            end
            
            % switch off these options
            options.onlykogeneassoc=false;
            options.deletegenes=false;
            options.ignorelistgenes={};
            options.dontkoess=false;
        end
    end
end


% check if reactions and genes to be ignored are present
if ~isempty(options.ignorelistrxns)
    rxnInModel=ismember(options.ignorelistrxns,model.rxns);
    if sum(rxnInModel)~=length(rxnInModel)
        % a banned reaction is not in the model
        warning(['Following reaction(s) to be ignored are not in model: ', char(join(options.ignorelistrxns(~rxnInModel),', ')) ])
    end
end
if ~isempty(options.ignorelistgenes)
    geneInModel=ismember(options.ignorelistgenes,model.genes);
    if sum(geneInModel)~=length(geneInModel)
        % a banned reaction is not in the model
        warning(['Following gene(s) to be ignored are not in model: ', char(join(options.ignorelistgenes(~geneInModel),', ')) ])
    end
end

% replace slashes, spaces and pluses in reaction/gene names since these characters
% will mess things up later
if options.deletegenes
    model.genes=strrep(model.genes,' ','');
    model.genes=strrep(model.genes,'/','\');
    model.genes=strrep(model.genes,'+','plus');
    if ~isempty(options.ignorelistgenes)
        options.ignorelistgenes=strrep(options.ignorelistgenes,' ','');
        options.ignorelistgenes=strrep(options.ignorelistgenes,'/','\');
    end
else
    model.rxns=strrep(model.rxns,' ','');
    model.rxns=strrep(model.rxns,'/','\');
    if ~isempty(options.ignorelistrxns)
        options.ignorelistrxns=strrep(options.ignorelistrxns,' ','');
        options.ignorelistrxns=strrep(options.ignorelistrxns,'/','\');
    end
    targetRxn=strrep(targetRxn,' ','');
    targetRxn=strrep(targetRxn,'/','\');
end

% check model is feasible and growth over threshold can occur
s=optimizeCbModel(model);
if s.stat~=1
    error('FBA failed on this model- model is probably infeasible')
elseif s.f==0
    error('Biomass reaction cannot carry flux')
elseif s.f<options.mingrowth
    error('Biomass reaction is below minimum flux threshold- please lower minGrowth')
end

% put model name in options structure so it is saved with parameters
if isfield(model,'modelID')
    options.modelname=model.modelID;
elseif isfield(model,'description')
    options.modelname=model.description;
end




%% CHECKING TARGET RXN

% check target reaction is in correct format
if isempty(targetRxn)
    error('No target reaction supplied')
end
if iscell(targetRxn)
    targetRxn=targetRxn{:};
end
if isstring(targetRxn)
    targetRxn=char(targetRxn);
end
if ~ischar(targetRxn)
    error('Target reaction is not a char')
end
if ~ismember(targetRxn,model.rxns)
    
    % check if the target relates to a metabolite in an exchange reaction
    excRxns=find(findExcRxns(model));
    [excMets,~]=find(model.S(:,excRxns)~=0);
    matchingMets=ismember(lower(model.mets(excMets)),lower(targetRxn));
    if isfield(model,'metNames')
        matchingMets=or(matchingMets, ismember( lower(model.metNames(excMets)), lower(targetRxn)) );
    end
    if isfield(model,'metFormulas')
        matchingMets=or(matchingMets, ismember(model.metFormulas(excMets), targetRxn));
    end
    
    if sum(matchingMets)==1
        % if target reaction refers to a metabolite in an exchange reaction, set exchange reaction as target
        [~,relevantExcRxns]=find(model.S(excMets(matchingMets),excRxns)~=0);
        if length(relevantExcRxns)==1
            targetRxn=model.rxns{excRxns(relevantExcRxns)};
        else
            disp('Target metabolite could refer to one of the following exchange reactions:')
            surfNet(model,model.rxns(excRxns(relevantExcRxns)))
            error('Unclear which exchange reaction should be optimised for the supplied metabolite- please enter a reaction instead of this metabolite')
        end
    else
        % target was not a metabolite in an exchange reaction- find if it is a metabolite and supply names of reactions that relate to this metabolite
        matchingMets=ismember(lower(model.mets),lower(targetRxn));
        if isfield(model,'metNames')
            matchingMets=or(matchingMets, ismember( lower(model.metNames), lower(targetRxn)) );
        end
        if isfield(model,'metFormulas')
            matchingMets=or(matchingMets, ismember(model.metFormulas, targetRxn));
        end
        if sum(matchingMets)==1
            disp(['Multiple reactions correspond to the metabolite: ', targetRxn])
            [~,matchingRxns]=find(model.S(matchingMets,:)~=0);
            surfNet(model,model.rxns(matchingRxns))
            error('Unclear which reaction should be optimised- please set the target as one of these reactions')
        elseif sum(matchingMets)>1
            disp('Unclear which of the following metabolites is the target:')
            disp(model.mets(matchingMets)')
            if isfield(model,'metNames')
                disp(model.metNames(matchingMets)')
            end
            if isfield(model,'metFormulas')
                disp(model.metFormulas(matchingMets)')
            end
            disp('Showing reactions that relate to the identified metabolites:')
            [~,matchingRxns]=find(model.S(matchingMets,:)~=0);
            surfNet(model,model.rxns(matchingRxns))
            error('Unclear which metabolite the target is- please set the target as a reaction that relates to one of these metabolites')
        else
            error(['Could not find anything in the model that corresponds to: ', targetRxn])
        end
    end
end

% check target reaction can carry flux
m=changeObjective(model,targetRxn);
s=optimizeCbModel(m);
if s.f==0 || s.stat~=1
    error('Target reaction cannot be active with current constraints')
elseif s.f<options.minprod
    error('Target reaction cannot exceed minimum flux threshold- consider adjusting constraints or minprod')
end
s=optimizeCbModel(m,'min');
if s.f<-options.tol
    warning('Reverse direction of target reaction is feasible- make sure that target reaction is written such that the forward reaction is the one that you want to happen. If it is not, multiply the column of the stoichiometric matrix (model.S) that corresponds to the target reaction by -1, save the updated model, and then start the algorithm again.')
    surfNet(model,targetRxn)
end

% put target reaction in options so it is saved with rest of parameters
options.targetreaction=targetRxn;

end

function output=checkInputType(input,type,name,inputSize)
% a function to check if a supplied input is in the right format and is the
% correct size


% check if no value supplied
if isempty(input)
    output=[];
    return
end

% if input is a cell and shouldn't be, and cell only contains one thing, set input to be the contents of the cell
if iscell(input) && ~strcmp(type,'cell')
    if isscalar(input)
        input=input{:};
    end
end

% check if variable is the correct type
switch type
    
    case 'numeric'
        if isnumeric(input)
            output=input;
        elseif ischar(input) || isstring(input)
            input=str2double(input);
            if isnan(input)
                output=[];
                warning(['options.',name, ' is not a number- default value will be used'])
                return
            else
                output=input;
            end
        else
            warning(['options.',name, ' is not a number- default value will be used'])
            output=[];
            return
        end
        
    case 'char'
        if ischar(input)
            output=input;
        else
            if isstring(input)
                output=char(input);
            else
                warning(['options.',name, ' is not a char or string- default value will be used'])
                output=[];
                return
            end
        end
        
    case 'cell'
        if iscell(input)
            output=input;
        else
            warning(['options.',name, ' is not a cell- default value will be used'])
            output=[];
            return
        end
        
    case 'double'
        if isdouble(input)
            output=input;
        elseif isnumeric(input) || islogical(input)
            output=double(input);
        else
            warning(['options.',name, ' is not a double- default value will be used'])
            output=[];
            return
        end
        
    case 'logical'
        if islogical(input)
            output=input;
        elseif isnumeric(input)
            if input==1
                output=true;
            elseif input==0
                output=false;
            else
                warning(['options.',name, ' is not a logical value- default value will be used'])
                output=[];
                return
            end
        elseif ischar(input) || isstring(input)
            input=str2double(input);
            if input==1
                output=true;
            elseif input==0
                output=false;
            else
                warning(['options.',name, ' is not a logical value- default value will be used'])
                output=[];
                return
            end
        else
            warning(['options.',name, ' is not a logical value- default value will be used'])
            output=[];
            return
        end
        
    case 'posInteger'
        if isnumeric(input)
            if floor(input)~=input
                output=[];
                warning(['options.',name, ' is not an integer- default value will be used'])
                return
            end
            if input<0
                output=[];
                warning(['options.',name, ' cannot be less than zero- default value will be used'])
                return
            end
            output=input;
        elseif ischar(input) || isstring(input)
            input=str2double(input);
            if isnan(input)
                output=[];
                warning(['options.',name, ' is not a number- default value will be used'])
                return
            else
                if floor(input)~=input
                    output=[];
                    warning(['options.',name, ' is not an integer- default value will be used'])
                    return
                end
                if input<0
                    output=[];
                    warning(['options.',name, ' cannot be less than zero- default value will be used'])
                    return
                end
                output=input;
            end
        else
            warning(['options.',name, ' is not a number- default value will be used'])
            output=[];
            return
        end
        
    case 'posNumeric'
        if isnumeric(input)
            output=input;
        elseif ischar(input) || isstring(input)
            input=str2double(input);
            if isnan(input)
                output=[];
                warning(['options.',name, ' is not a number- default value will be used'])
                return
            else
                output=input;
            end
        else
            warning(['options.',name, ' is not a number- default value will be used'])
            output=[];
            return
        end
        if output<0
            warning(['options.',name, ' cannot be less than 0- default value will be used'])
        end
        
end

% check if variable is the correct size
if sum(size(output)>inputSize)~=0
    if sum(size(output)>inputSize([2,1]))==0
        % input is an appropriate size but is in wrong orientation
        output=output';
    else
        % input is too big- default will be used in case an inappropriately sized input causes problems later
        warning(['options.',name, ' is too large- it should be a maximum of: ',num2str(inputSize), ' and was actually: ' num2str(size(output)), ' default value will be used instead'])
        output=[];
    end
end


end