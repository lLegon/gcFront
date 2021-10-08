function [FBAsolution RobustKnocksolution] = myFBA(model, p)
%function [FBAsolution RobustKnocksolution] = myFBA(model, p)
%Modified Flux-Balance-Analysis for maximizing objectivem and maximizing/minimizing outer objective function 
%Written by Andras Hartmann <andras.hartmann@gmail.com> based on COBRA optimizeCBModel function
%
% First solves LP problem of the form: 
%                             f = max   model.c'*v
%                                 s.t.  S*v = b         : y
%                                       lb <= v <= ub   : w
%
%
% Then solves th LP problem:
%                                 min   p'v
%                                 s.t.  S*v = b
%                                       c'v = f
%                                       lb <= v <= ub
%                           
%
%INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%
% p  the desired product objective
%
%OUTPUT
% FBAsolution and RobustKnocksolution
%   f         Objective value
%   x         Primal
%   y         Dual
%   w         Reduced costs
%   s         Slacks
%   stat      Solver status in standardized form
%              1   Optimal solution
%              2   Unbounded solution
%              0   Infeasible
%             -1  No solution reported (timelimit, numerical problem etc)
%   origStat  Original status returned by the specific solver


gurobi_error = 1;

while gurobi_error

try

gurobi_error = 0;
primalOnlyFlag = 0;
printLevel = 0;
[nMets,nRxns] = size(model.S);

%% Setting up the first LP problem
LPproblem.A = model.S;
LPproblem.c = model.c;
LPproblem.b = zeros(size(model.S,1),1);
LPproblem.lb = model.lb;
LPproblem.ub = model.ub;
LPproblem.osense = -1; %maximizing
LPproblem.csense(1:nMets,1) = 'E';

%%Double check that all inputs are valid:
%if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
%    warning('invalid problem');
%    return;
%end
%should be validating, but faster if not

%t1 = clock;
%also no clock for faster solution

%% Solve initial LP
    solution = solveCobraLP(LPproblem);

if (solution.stat ~= 1)||(solution.obj == 0)
    % check if initial solution was successful.
    %{
    if printLevel>0
        warning('Optimal solution was not found');
    end
    %}
    FBAsolution.f = 0;
    FBAsolution.x = [];
    FBAsolution.stat = solution.stat;
    FBAsolution.origStat = solution.origStat;
    FBAsolution.solver = solution.solver;
    %FBAsolution.time = etime(clock, t1);
    
    RobustKnocksolution = FBAsolution;

    return;
end

objective = solution.obj; % save for later use.

%% Store results
if (solution.stat == 1)
    %solution found.
    FBAsolution.x = solution.full(1:nRxns);
    
    if isfield(solution,'dual')
        if ~isempty(solution.dual)
            %dont include dual variable to additional constraint
            solution.dual=solution.dual(1:end-1,1);
        end
    end
    
    %this line IS necessary.
    FBAsolution.f = model.c'*solution.full(1:nRxns); %objective from original optimization problem.
    if abs(FBAsolution.f - objective) > .01
            error('myFBA.m: minimizing Euclidean norm did not work')
    end
    
    %if (~primalOnlyFlag && allowLoops && any(~minNorm)) % LP rcost/dual only correct if not doing minNorm
    % LP rcost/dual are still meaninful if doing, one simply has to be aware that there is a
    % perturbation to them the magnitude of which depends on norm(minNorm) - Ronan   
    if (~primalOnlyFlag)
        FBAsolution.y = solution.dual;
        FBAsolution.w = solution.rcost;
    end
else
    %some sort of error occured.
    if printLevel>0
        warning('Optimal solution was not found');
    end
    FBAsolution.f = 0;
    FBAsolution.x = [];
end

FBAsolution.stat = solution.stat;
FBAsolution.origStat = solution.origStat;
FBAsolution.solver = solution.solver;
%FBAsolution.time = etime(clock, t1);

%%setting up second LP
LPproblem2=LPproblem;
LPproblem2.A = [LPproblem.A; LPproblem.c'];%new constraint a row with the old objective
LPproblem2.b = [LPproblem.b; solution.full(LPproblem.c~=0)];
LPproblem2.csense(end+1) = 'E';
LPproblem2.c = p; %minimize the target production
LPproblem2.osense = 1; %I said minimize

%% Solve initial LP
    solution = solveCobraLP(LPproblem2);

if (solution.stat ~= 1) % check if initial solution was successful.
    %{
    if printLevel>0
        warning('Optimal solution was not found');
    end
    %}
    RobustKnocksolution.f = 0;
    RobustKnocksolution.x = [];
    RobustKnocksolution.stat = solution.stat;
    RobustKnocksolution.origStat = solution.origStat;
    RobustKnocksolution.solver = solution.solver;
    %RobustKnocksolution.time = etime(clock, t1);
    return;
end

objective = solution.obj; % save for later use.

%% Store results
if (solution.stat == 1)
    %solution found.
    RobustKnocksolution.x = solution.full(1:nRxns);
    
    if isfield(solution,'dual')
        if ~isempty(solution.dual)
            %dont include dual variable to additional constraint
            solution.dual=solution.dual(1:end-1,1);
        end
    end
    
    %this line IS necessary.
    RobustKnocksolution.f = p'*solution.full(1:nRxns); %objective from original optimization problem.
    if abs(RobustKnocksolution.f - objective) > .01
            error('myFBA.m: minimizing Euclidean norm did not work')
    end
    
    %if (~primalOnlyFlag && allowLoops && any(~minNorm)) % LP rcost/dual only correct if not doing minNorm
    % LP rcost/dual are still meaninful if doing, one simply has to be aware that there is a
    % perturbation to them the magnitude of which depends on norm(minNorm) - Ronan   
    if (~primalOnlyFlag)
        RobustKnocksolution.y = solution.dual;
        RobustKnocksolution.w = solution.rcost;
    end
else
    %some sort of error occured.
    if printLevel>0
        warning('Optimal solution was not found');
    end
    RobustKnocksolution.f = 0;
    RobustKnocksolution.x = [];
end

RobustKnocksolution.stat = solution.stat;
RobustKnocksolution.origStat = solution.origStat;
RobustKnocksolution.solver = solution.solver;
%RobustKnocksolution.time = etime(clock, t1);

catch 
    disp('ooops gurobi experimented some mailfunction')
    gurobi_error = 1;
end
end
