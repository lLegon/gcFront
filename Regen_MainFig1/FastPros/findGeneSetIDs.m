function geneSetID = findGeneSetIDs(model,geneList)
%findGeneSetIDs Find gene set numbers in a model
%
% geneID = findGeneSetIDs(model,geneList)
%
%INPUTS
% model     Reduced COBRA model structure created by reduceModelForFP
% geneList  List of genes
%
%OUTPUT
% geneSetID    List of gene IDs corresponding to geneList
%
% 
% Aug. 5th, 2013    Satoshi OHNO

if (iscell(geneList))
    [~,geneSetID] = ismember(geneList,model.geneSets);    
else
    geneSetID = find(strcmp(model.geneSets,geneList));
    if (isempty(geneSetID))
        geneSetID = 0;
    end
    if (length(geneSetID) > 1)
        geneSetID = geneSetID(1);
    end
end