function model_ = delReaction(thismodel, rxnNo)
%simply delete a reaction
    model_ = thismodel;
    model_.lb(rxnNo)=0;
    model_.ub(rxnNo)=0;
end
