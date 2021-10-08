function results_matched = matchCompr2OrigMod(results,model_s,model,probOpts)
% Match found reactions of compressed model to original model

reacidx     = model_s.reacidx;
KORxnNum    = results.KORxnNum;
KONumMatch  = {};

%% Differentiante multiple equal solutions (lumped reactions)
numIntStrat     = 1;
KONumMatch      = cell(2,size(KORxnNum,1));
for i=1:size(KORxnNum,1)
    KONumMatch{1,i}     = find(reacidx(:,KORxnNum(i))); 
    numRxns             = length(KONumMatch{1,i});
    numIntStrat         = numIntStrat*numRxns;
    KONumMatch{2,i}     = numRxns;
end

KoNumStore  = zeros(size(KONumMatch,2),numIntStrat);
for i=1:size(KONumMatch,2)
    numRxns     = KONumMatch{2,i};
    rxnIdxWrite = KONumMatch{1,i};
    x           = numIntStrat/numRxns;
    for j=1:x
        initIdx     = ((j-1)*numRxns)+1;
        KoNumStore(i,initIdx:(initIdx+numRxns-1))  = rxnIdxWrite';
    end   
end


%% Rewrite results file
% KOs
results_matched.KORxnNum            = KoNumStore;
results_matched.ClusteredKORxns     = KONumMatch(1,:);

for i=1:size(KoNumStore,2)
    results_matched.KORxnNames{i}   = model.rxnNames(KoNumStore(:,i));
    KOs(:,i)                        = zeros(size(model.rxns,1),1);
    KOs(KoNumStore(:,i),i)            = 1;
end
results_matched.KOs             = KOs;


% Fluxes
results_matched.fluxes          = reacidx*results.fluxes;
if isfield(results,'FVAfluxes')
    results_matched.FVAfluxes       = reacidx*results.FVAfluxes;
end

end