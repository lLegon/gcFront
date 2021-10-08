function solution = optimizeCbModel2(LPproblem)
% ------------------------------------------------------------------------
% Function modified by Laurence Legon on 30/09/2021
% 
% This code is a copy of the second half half of the COBRA toolbox function "optimizeCbModel"
% Accepts the LP problem from "convertModel" and solves it. Can be used to speed up code by
% avoiding redundant conversion of models into LP problem format
%
% .. Authors of the original optimizeCbModel function:
%       - Markus Herrgard       9/16/03
%       - Ronan Fleming         4/25/09  Option to minimises the Euclidean Norm of internal
%                                        fluxes using 'cplex_direct' solver
%       - Ronan Fleming         7/27/09  Return an error if any imputs are NaN
%       - Ronan Fleming         10/24/09 Fixed 'E' for all equality constraints
%       - Jan Schellenberger             MILP option to remove flux around loops
%       - Ronan Fleming         12/07/09 Reworked minNorm parameter option to allow
%                                        the full range of approaches for getting
%                                        rid of net flux around loops.
%       - Jan Schellenberger    2/3/09   fixed bug with .f being set incorrectly
%                                        when minNorm was set.
%       - Nathan Lewis          12/2/10  Modified code to allow for inequality
%                                        constraints.
%       - Ronan Fleming         12/03/10 Minor changes to the internal handling of
%                                        global parameters.
%       - Ronan Fleming         14/09/11 Fixed bug in minNorm with negative
%                                        coefficient in objective
%       - Minh Le               11/02/16 Option to minimise the cardinality of
%                                        fluxes vector
%       - Stefania Magnusdottir 06/02/17 Replace LPproblem2 upper bound 10000 with Inf
%       - Ronan Fleming         13/06/17 Support for coupling C*v<=d




%%

[printLevel,primalOnlyFlag] = getCobraSolverParams('LP',{'printLevel','primalOnly'});

% size of the stoichiometric matrix
[nMets,nRxns] = size(LPproblem.A);


if isfield(LPproblem,'C')
    nCtrs = size(LPproblem.C,1);
end

if isfield(LPproblem,'E')
    nVars = size(LPproblem.E,2);
end

t1 = clock;



% Solve initial LP

    solution = solveCobraLP(LPproblem);


global CBT_LP_SOLVER
if strcmp(CBT_LP_SOLVER,'mps')
    solution=solution;
    return;
else
    if (solution.stat ~= 1) % check if initial solution was successful.
        if printLevel>0
            warning('Optimal solution was not found');
        end
        solution.f = 0;
        solution.x = [];
        solution.stat = solution.stat;
        solution.origStat = solution.origStat;
        solution.solver = solution.solver;
        solution.time = etime(clock, t1);
        return;
    end
end

objective = solution.obj; % save for later use.

[nTotalConstraints,nTotalVars] = size(LPproblem.A);


% Store results
if (solution.stat == 1)
    % solution found. Set corresponding values
    solution.x = solution.full(1:nRxns);
    solution.v = solution.x;
    % handle the objective, otherwise another objective value could be 
    % returned and we only want to return the value of the defined
    % model objective
    if isfield(LPproblem,'E')
        solution.vars_v = solution.full(nRxns+1:nRxns+nVars);        
        solution.f = LPproblem.c'*solution.v + LPproblem.evarc' * solution.vars_v; % We need to consider the 
    else
        solution.f = LPproblem.c'*solution.full(1:nRxns); %objective from original optimization problem.
    end
    % Check objective quality
    if abs(solution.f - objective) > .01
        error('optimizeCbModel.m: minimizing Euclidean norm did not work')
    end
    
    % handle the duals, reducing them to fields in the model.
    if isfield(solution,'dual')
        if ~isempty(solution.dual)
            if isfield(LPproblem,'C')
                solution.ctrs_y = solution.dual(nMets+1:nMets+nCtrs,1);
            end
            solution.dual=solution.dual(1:nMets,1);            
        end    
    end            
    
    % handle reduced costs 
    if isfield(solution,'rcost')
        if ~isempty(solution.rcost)
            if isfield(LPproblem,'E')
                solution.vars_w = solution.rcost(nRxns+1:nRxns+nVars,1);
            end
            solution.rcost=solution.rcost(1:nRxns,1);            
        end
    end     
    
    % handle slacks
    if isfield(solution,'slack')
        if ~isempty(solution.slack)
            if isfield(LPproblem,'C')
                solution.ctrs_s = solution.slack(nMets+1:nMets+nCtrs,1);
            end
            solution.slack=solution.slack(1:nMets,1);            
        end
    end     
    
    %if (~primalOnlyFlag && allowLoops && any(~minNorm)) % LP rcost/dual only correct if not doing minNorm
    % LP rcost/dual are still meaninful if doing, one simply has to be aware that there is a
    % perturbation to them the magnitude of which depends on norm(minNorm) - Ronan
    if (~primalOnlyFlag)
        solution.y = solution.dual;          
        solution.w = solution.rcost; 
        solution.s = solution.slack;
    end
    fieldOrder = {'full';'obj';'rcost';'dual';'slack';'solver';'algorithm';'stat';'origStat';'time';'basis';'vars_v';'vars_w';'ctrs_y';'ctrs_s';'f';'x';'v';'w';'y';'s'};
    % reorder fields for better readability
    currentfields = fieldnames(solution);
    presentfields = ismember(fieldOrder,currentfields);
    absentfields = ~ismember(currentfields,fieldOrder);
    solution = orderfields(solution,[currentfields(absentfields);fieldOrder(presentfields)]);
else
    %some sort of error occured.
    if printLevel>0
        warning('Optimal solution was not found');
    end
    solution.f = 0;
    solution.x = [];
    solution.v = solution.x;
end
solution.time = etime(clock, t1);







