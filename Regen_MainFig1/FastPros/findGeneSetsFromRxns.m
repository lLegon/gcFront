function [geneSetList]=findGeneSetsFromRxns(model,reactions)
%findGeneSetsFromRxns make a gene set list of the genes that correspond to the 
%selected reactions
%
% [geneList]=findGeneSetsFromRxns(model,reactions)
% 
%INPUTS
% model         Reduced COBRA model structure created by reduceModelForFP
% reactions     Reactions to find the corresponding genes for
%
%OUTPUT
% geneList      List of genes corresponding to reactions
% 
% Aug. 5th, 2013    Satoshi OHNO

rxnInd = findRxnIDs(model, reactions);
RxnNotInModel = find(rxnInd==0);
if ~isempty(RxnNotInModel)
    for i = 1:length(RxnNotInModel)
        display(cat(2,'The reaction "', reactions{RxnNotInModel(i)},'" is not in your model!'));
        
    end
end
rxnInd(RxnNotInModel) = [];
reactions(RxnNotInModel) = [];
for i = 1:length(rxnInd)
    geneSetList(i) = model.geneSets(find(model.geneSetRxnMat(:,rxnInd(i))));
end

geneSetList = columnVector(geneSetList);