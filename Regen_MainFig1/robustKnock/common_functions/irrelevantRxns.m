function [irrel relevantRxns relevantRxnIDs] = irrelevantRxns2(model,biomassRxn,userRxns)
%function [irrel relevantRxns relevantRxnIDs] = irrelevantRxns(model,biomassRxn)
%
%This function is for pre-processing. It determines a list of "irrelevant" reactions that should not be considered for KO.
% Irrelevant reactions are the essential, transport and blocked blocked reactions. 

%obsoleted, should not be called with one argument
%if nargin < 2
%    biomassRxn = 'biomass_a';
%end

%% Initialization

%Just in case objective would be different
cobra_model=changeObjective(model,biomassRxn);


%Find biomass rxn
biomassrxn = find(ismember(model.rxns,biomassRxn)); %301
% lower flux limit for reaction to be considered lethal 
fluxTol = 1e-9; 

%% Lethal Rxns
[~,~,~,~,~,Flux_FBA] = singleRxnDeletion(model,'FBA');
irrel.lethalRxns = Flux_FBA( biomassrxn,: )< fluxTol; 
%irrel.ind = find (Flux_FBA( biomassrxn,: )< fluxTol ); 

%% Transport and synthetic reactions
%sink
irrel.tr_synthRxns = (cellfun(@isempty, ...
        strfind(model.rxns, 'sink_'))'==0);
%export
irrel.tr_synthRxns = (irrel.tr_synthRxns|(cellfun(@isempty, ...
        strfind(model.rxns, '(e)'))'==0));
irrel.tr_synthRxns = (irrel.tr_synthRxns|(cellfun(@isempty, ...
        strfind(model.rxns, 'export'))'==0));
irrel.tr_synthRxns = (irrel.tr_synthRxns|(cellfun(@isempty, ...
    strfind(model.rxns, 'EX_'))'==0));
%TODO: biomass (is this ok this way?)
irrel.tr_synthRxns = (irrel.tr_synthRxns|(cellfun(@isempty, ...
        strfind(model.rxns, 'ass'))'==0));
    
%% Blocked reactions
irrel.blockedRxns = ismember(model.rxns, findBlockedReaction(model))';

irrelevant = irrel.blockedRxns|irrel.tr_synthRxns|irrel.lethalRxns;

irrel.ind=find(irrelevant)
irrel.names=model.rxns(irrel.ind);

% ADDED CODE SO REACTIONS NOT ON USER CONSIDERED LIST WILL BE IGNORED


relevantRxnIDs = find(~irrelevant);
relevantRxns = model.rxns(relevantRxnIDs);

%% Save
%save('irrelevantRxns','irrel');
end

