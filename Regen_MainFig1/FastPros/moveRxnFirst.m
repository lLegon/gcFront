function [model] = moveRxnFirst(model,movespot)
% The function moves a reaction from one spot in the network to first.
%
% [model] = moveRxnFirst(model,movespot)
%
%INPUTS
% model         COBRA model structure
% movespot     The reaction number to move to first
% 
%OUTPUTS
% model         COBRA toolbox model structure with moved reaction
% 
% Aug. 5th, 2013    Satoshi OHNO

oldmodel = model;
restspot = setdiff(1:length(model.rxns),movespot);
movespot = columnVector(movespot)';

model.lb = oldmodel.lb([movespot,restspot]);
model.ub = oldmodel.ub([movespot,restspot]);
model.c = oldmodel.c([movespot,restspot]);

if isfield(model,'rxns')
    model.rxns = oldmodel.rxns([movespot,restspot]);
end

if isfield(model,'rxnNames')
    model.rxnNames = oldmodel.rxnNames([movespot,restspot]);
end

if isfield(model,'subSystems')
    model.subSystems = oldmodel.subSystems([movespot,restspot]);
end

if isfield(model,'rules')
    model.rules = oldmodel.rules([movespot,restspot]);
end

if isfield(model,'grRules')
    model.grRules = oldmodel.grRules([movespot,restspot]);
end

if isfield(model,'rev')
    model.rev = oldmodel.rev([movespot,restspot]);
end

model.S = oldmodel.S(:,[movespot,restspot]);

model.rxnGeneMat = oldmodel.rxnGeneMat([movespot,restspot],:);
