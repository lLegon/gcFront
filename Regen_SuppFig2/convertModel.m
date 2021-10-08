function LPproblem=convertModel(model, osenseStr, minNorm, allowLoops, zeroNormApprox)
% Function modified by Laurence Legon (30/09/21)
% Essentially, this is the first half of the COBRA toolbox function "optimizeCbModel"
% Converts a COBRA model into a LP problem that can be solved by "optimizeCbModel2". 
% Can be used to speed up code by avoiding redundant conversion of models into LP problem format
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
if exist('osenseStr', 'var') % Process arguments and set up problem
    if isempty(osenseStr)
        model.osenseStr = 'max';
    else
        model.osenseStr = osenseStr;
    end
else
    if isfield(model, 'osenseStr')
        model.osenseStr = model.osenseStr;
    else
        model.osenseStr = 'max';
    end
end


%use global solver parameter for printLevel


LPproblem = buildLPproblemFromModel(model);

%Double check that all inputs are valid:
if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
    warning('invalid problem');
    return;
end


end