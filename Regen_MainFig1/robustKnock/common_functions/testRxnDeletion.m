function [delRxn,fluxSolution] = testRxnDeletion(model,evaluate,numDel,rxnList, timeLimit, startTime, endTol)
% Modified by Laurence Legon (30/09/21) to allow termination after a GC
% design is found, or after a time limit has been exceeded
%
%function [delRxn,fluxSolution] = testRxnDeletion(model,evaluate,numDel,rxnList)

if (nargin < 4)
    error('testRxnDeletion needs at least 3 arguments')
end

if (nargin < 4)
    rxnList = 1:(size(model.rxns));
else
    if (isempty(rxnList))
        rxnList = 1:(size(model.rxns));
    end
end

tol = 10^-9;
if toc(startTime)>timeLimit
    delRxn=[];
    fluxSolution=[];
    return
end

nDelRxns = length(rxnList);
if numDel == 1
    delRxn = zeros(nDelRxns,1);
    fluxSolution = zeros(nDelRxns,2);
    parfor i = 1:nDelRxns
        changeCobraSolver('gurobi','all',0);
        delRxn(i) = rxnList(i);
        fluxSolution(i,:) = evaluate(delReaction(model, rxnList(i)));
    end;
    return;
else
    delRxn = [];
    fluxSolution = [];
    for i = 1:nDelRxns-numDel+1
        if ~isempty(fluxSolution) && max(fluxSolution(:,2))>endTol
            % SKIP IF ENDTOL HAS BEEN EXCEEDED
            continue
        end
        %Uncomment this if you would like to have some status
        %if numDel == 3
        %    i
        %end
        [delRxn_,fluxSolution_] = testRxnDeletion(delReaction(model, rxnList(i)), evaluate, numDel-1, rxnList(i+1:end), timeLimit, startTime, endTol);
        
        if isempty(delRxn_) % IF statement added so the empty outputs created after timelimit exceeded aren't read
            continue
        end
        
        delRxn_=[ones(size(delRxn_,1),1)*rxnList(i) delRxn_];
        %Only those solutions where biomass larger than tolerance
        
        delRxn = [delRxn; delRxn_( find(fluxSolution_(:,2)>tol), :) ];
        fluxSolution = [fluxSolution; fluxSolution_( find(fluxSolution_(:,2)>tol), :)];
        
        
        
    end;
end;

end

