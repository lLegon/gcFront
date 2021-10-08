function rxnID = findCombinedRxnIDs(model,rxnList)
%findCombined Find combined reaction numbers in a reduced model
%
% rxnID = findCombinedRxnIDs(model,rxnList)
%
%INPUTS
% model     Reduced COBRA model strcture created by reduceModelForFP
% rxnList   List of reactions in original model
%
%OUTPUT
% rxnID     IDs for reactions corresponding to rxnList in the reduced model
%
% Aug. 5th, 2013    Satoshi OHNO


rxnID = zeros(length(rxnList),1);
for i = 1: length(rxnList)
    [temp,~]= find(strcmp(rxnList(i), model.rxnAssociations));
    if ~isempty(temp)
        rxnID(i) = temp;
    end
end


