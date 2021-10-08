function [rxnSetList]=findRxnSetsFromGeneSetIDs(model,geneSetIDs)
%findRxnSetsFromGeneSetIDs make a reaction set list of the gene sets that correspond to the 
%selected gene set numbers
%
% [rxnSetList]=findRxnSetsFromGeneSetIDs(model,geneSetIDs)
% 
%INPUTS
% model         Reduced COBRA model structure created by reduceModelForFP
% geneSetIDs    IDs of gene sets to find the corresponding rxn sets for
%
%OUTPUT
% rxnSetList      List of rxn sets corresponding to reactions
% 
% Aug. 5th, 2013    Satoshi OHNO


rxnSetList = cell(size(geneSetIDs));
for i = 1 : numel(geneSetIDs)
    rxnSetList(i) = model.geneSetAssocRxns(geneSetIDs(i));
end
