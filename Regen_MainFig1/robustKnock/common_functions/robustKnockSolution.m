function fluxSol = robustKnockSolution(model,p)
%function fluxSol = robustKnockSolution(model,p)
%[~, rb] = myFBA(model, p);
[fb, rb] = myFBA(model, p);
if (~isempty(rb.x))
    fluxSol = [rb.x(find(model.c)) rb.x(find(p))];
    return;
end
fluxSol = zeros(1,2);
