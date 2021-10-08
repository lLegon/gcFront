function [irrel relevantRxns] = irrelevantRxns(model,biomassRxn)
%function irrel = irrelevantRxns(model)
%This function creates a mat vector with the indexes of all "irrelevant"
%reactions: essential, transport and blocked. So the ones we don't want to
%KO

if nargin < 2
    biomassRxn = 'biomass_LLA';
end

%% Init
%Find biomass rxn
biomassrxn = find(ismember(model.rxns,biomassRxn)); %301
% lower flux for reaction to be lethal 
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
%biomass (is this ok this way?)
irrel.tr_synthRxns = (irrel.tr_synthRxns|(cellfun(@isempty, ...
        strfind(model.rxns, 'ass'))'==0));
    
%% Blocked reactions
irrel.blockedRxns = ismember(model.rxns, findBlockedReaction(model))';

irrelevant = irrel.blockedRxns|irrel.tr_synthRxns|irrel.lethalRxns;

irrel.ind=find(irrelevant)
irrel.names=model.rxns(irrel.ind);

relevantRxns = model.rxns(find(~irrelevant));

%% Save
save('irrelevantRxns_Lactis','irrel');
end

