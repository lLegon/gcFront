function [eRxns,resTest] = essentialRxns(S_s,lb_s,ub_s,cBm,B)
% Determine essential reactions or combinations of reactions which KO leads
% to "dead" cells

numRxns_s       = size(S_s,2);
numRxns         = size(B,1);
numConstr       = size(S_s,1);

I           = eye(numRxns_s);

%% Setup general LP problem
Av      = [S_s;-S_s;I;-I];

b       = [zeros(2*numConstr,1);ub_s;-lb_s];

% cBm     = [cBm;zeros(numRxns,1)];

%% Test model
% Setup gurobi model
gurMod.A        = sparse(Av);
gurMod.obj      = cBm;
gurMod.sense    = '<';
gurMod.rhs      = b;
gurMod.lb       = lb_s;
gurMod.ub       = ub_s;

vtype                       = cell(numRxns_s,1);
vtype(:)                    = {'C'};
gurMod.vtype                = char(vtype);

gurMod.modelsense           = 'max';

% Parameter
params.FeasibilityTol   = 1e-05;
params.Presolve         = 1;        % Dissable Presolve (Disturbes Solution)
params.Threads          = 1;
params.OutputFlag       = 0;

resTest     = gurobi(gurMod,params);


%% Loop through reactions
params.OutputFlag       = 0;

B       = logical(B);
eRxns   = zeros(numRxns,1);
gurModStore     = cell(numRxns,1);
% Prepare models
for r=1:numRxns
    % Set KOs
    lb_parse    = lb_s;
    ub_parse    = ub_s;
    lb_parse(B(r,:))    = 0;
    ub_parse(B(r,:))    = 0;
    
    % Update model
    gurMod.lb       = lb_parse;
    gurMod.ub       = ub_parse;
    
    gurModStore{r}  = gurMod;     
end

parfor r=1:numRxns
    % Attempt solving LP problem
    res     = gurobi(gurModStore{r},params);
    
    if strcmp(res.status,'INFEASIBLE')
        eRxns(r)    = 1;
    end

end

end