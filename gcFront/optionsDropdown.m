function [dd,event]=optionsDropdown(dd,lbl,optionFig,textEntry)
% code that responds to the user selecting something in the options entry
% dropdown menu


% display information on what the variable picked is, and change the
% variable that the textbox links to
switch dd.Value
    
    case 'solver'
        lbl.Text={'Name of the linear programming solver to be used. Leave blank to use', 'currently initialised solver. If no solver initialised, will use COBRA', 'toolbox default solver.'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'solver');
        
    case 'biomassrxn'
        lbl.Text={'Name of the biomass reaction. Leave blank to use the current', 'model objective.'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'biomassrxn');
        
    case 'tol'
        lbl.Text={'Tolerance to mathematical error. Fluxes below this value will be', 'ignored.'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'tol');
        
    case 'tiltval'
        lbl.Text={'The coefficient used for objective value tilting (so minimum product', 'flux is found).'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'tiltval');
        
    case 'shiftval'
        lbl.Text={'The size of the reduction in growth rate that is used to calculate', 'shadow price'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'shiftval');
        
    case 'skipreduction'
        lbl.Text={'Enter 0 to reduce size of model by removing inactive reactions and', 'pooling linear pathways, or 1 to skip this'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'skipreduction');
        
    case 'mingrowth'
        lbl.Text={'Ignore designs with growth below this value'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'mingrowth');
        
    case 'minprod'
        lbl.Text={'Do not consider deletions if they will lower product synthesis below', 'this threshold'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'minprod');
        
    case 'removeredundancy'
        lbl.Text={'Enter 1 to remove redundant reactions from identified designs, or 0', 'to skip this'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'removeredundancy');
        
    case 'saveresults'
        lbl.Text={'Enter 1 to save parameters and designs identified, or 0 to skip this'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'saveresults');
        
    case 'maxknockouts'
        lbl.Text={'Maximum number of knockouts that a design may have'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'maxknockouts');
        
    case 'deletegenes'
        lbl.Text={'Enter 1 to search for gene knockouts, or 0 to search for reaction', 'knockouts'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'deletegenes');
        
    case 'ignorelistrxns'
        lbl.Text={'Enter the names of any reactions that should not be knocked out, with', 'commas inbetween reaction names (e.g. RXN1, RXN2, RXN3)','Reaction names are case sensitive'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'ignorelistrxns');
        
    case 'ignorelistgenes'
        lbl.Text={'Enter the names of any genes that should not be knocked out, with', 'commas inbetween gene names (e.g. GENE1, GENE2, GENE3).','Gene names are case sensitive'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'ignorelistgenes');
        
    case 'dontkoess'
        lbl.Text={'Enter 1 to prevent reactions from being knocked out if they can only', 'be knocked out by knocking out essential genes'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'dontkoess');
        
    case 'onlykogeneassoc'
        lbl.Text={'Enter 1 to only knock out reactions if they are gene associated, or 0', 'if they are not'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'onlykogeneassoc');
        
    case 'mutationrate'
        lbl.Text={'Average number of changes that are made to each design during', 'the mutation step of the GA'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'mutationrate');
        
    case 'popsize'
        lbl.Text={'Size of the GA population'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'popsize');
        
    case 'genlimit'
        lbl.Text={'Maximum number of generations that the GA should run for'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'genlimit');
        
    case 'timelimit'
        lbl.Text={'Maximum duration of the GA (in seconds)'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'timelimit');
        
    case 'fitnesslimit'
        lbl.Text={'Terminate algorithm if a design is discovered with greater than this', 'product synthesis'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'fitnesslimit');
        
    case 'spreadchangelimit'
        lbl.Text={'Terminate algorithm if spread of designs is below this value for', 'stallgenlimit generations'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'spreadchangelimit');
        
    case 'stallgenlimit'
        lbl.Text={'Terminate algorithm if spread of designs is below spreadchangelimit', 'for this many generations'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'stallgenlimit');
        
    case 'plotinterval'
        lbl.Text={'Number of generations that should pass before the plot showing all', 'the designs should be updated'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'plotinterval');
        
    case 'newredundantremoval'
        lbl.Text={'If redundant KOs are to be removed, enter 1 to use new (faster)', 'removal of redundancy, or 0 to remove redundancy as was done in', 'the gcFront paper'};
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'newredundantremoval');
    
    case 'maxreductionsize'
        lbl.Text={'Designs with more KOs than this parameter will not have redundant', 'KOs fully removed, to reduce the number of combinations to analyse.', 'Note- only works with newredundantremoval method!' };
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'maxreductionsize');
        
    otherwise
        lbl.Text='Unrecognised option';
        textEntry.ValueChangedFcn=@(dd,event)updateOptionsFromFig(dd,event,'');
        textEntry.Value='';
        return
        
end
    

% populate the textbox with the current value
textEntry.Value=char(eval(['optionFig.UserData.',dd.Value]));




end