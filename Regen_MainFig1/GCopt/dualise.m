function [d_Av, d_Ay, d_b, d_c, lb_d, ub_d] = dualise(Av, Ay, b, c, solver, indLR, indLC, genOpts)
% Takes matrices of the primal LP problem and transform them into the dual
% system
%
% Primal System:
% min -cT v
% Av v +  Ay y -  b <= 0
%
% Output matrices:
% d_Av, d_Ay, d_b, d_S according to 
% d_Av d_v + d_Ay d_y <= d_b
%
%
%%
fprintf('--> Dualise LP-problem \n')

numIV       = size(Ay,2);       % Number of integer parameter (number of reactions before the split)(q)
numDVy      = size(Ay,1);       % Number of flux constraints (of dual variables)
numDVv      = size(Av,1);       % Number of flux constraints (of dual variables)
numPV       = size(Av,2);       % Number of primal variables (fluxes after the split)


Av_sprs     = sparse(Av);


%% Define numbers of auxiliary variables z
numAV       = size(indLR,1);


%% Do a matrices check
if numDVv ~= numDVy
    error('Row dimensions of Av and Ay do not match!')
end

numDV       = numDVy;       % Number of dual variables


%% Linearise bi-linear terms y*lambda
maxLFac = genOpts.avBnd;
maxL    = maxLFac*ones(numDV,1);

% Preallocate constraint matrices
Dv1_d       = sparse(numAV,numDV);
Dv2_d       = sparse(numAV,numDV);

Dy1_d       = sparse(numAV,numIV);
Dy2_d       = sparse(numAV,numIV);

b1_d        = ones(numAV,1);

c_0         = zeros(numAV,1);

% indLambda   = [indLambda1;indLambda2];

%%
% Gurobi Solver Options
gurModel.A           = -Av_sprs';         % linear constraint matrix
gurModel.lb          = zeros(numDV,1);   % lower bound vector
gurModel.sense       = '<';
% gurModel.rhs         = zeros(numPV,1);   % right-hand side vector linear constraints
gurModel.rhs         = full(-c);   % right-hand side vector linear constraints
gurModel.vtype       = 'C';              % variable type (continuous)
gurModel.modelsense  = 'max';

params.InfUnbdInfo   = 0;
params.OutputFlag    = 0;
params.Presolve      = 0;
params.Threads       = 1;

% Define options for each objective
gurModelStruct          = cell(numDV);
for r=1:numDV
    c_r                 = zeros(numDV,1);
    c_r(r,1)            = 1;
    gurModel.obj        = c_r;
    gurModelStruct{r}   = gurModel;
end

% % Write constraint matrices for auxiliary variables z
zCount      = 1;
solverLP    = solver.LP;
for r=1:numDV
    
    % Find upper bound for dual variable
%     c_r             = zeros(numDV,1);
%     c_r(r,1)        = 1;

    switch solverLP
        case 0      % Gurobi Solver
            
            res     = gurobi(gurModelStruct{r},params);

            % Check if solution is feasible
            if strcmp(res.status,'OPTIMAL')
%                 disp('Feasible solution')
                maxL(r,1)   = res.objval;
            elseif strcmp(res.status,'UNBOUNDED')
%                 disp(res.status)
            else
                error('Infeasible LP-problem during auxiliary variable maximisation!')
            end
    end
end
    
for r=1:numDV    
    [idx,~]         = find(indLR(:,1)==r,1);
    if ~isempty(idx)
        for z=1:length(idx)
            % Write constraint matrices for auxiliary variables z 
            dvIdx       = indLR(idx(z),1);       % Index of dual variable to be replaced
            zIdx        = indLC(idx(z),1);
            val         = Ay(dvIdx,zIdx);

            % Check val
            if val==0
                error('Assignment of auxiliary variables z failed!')
            end

            Dv1_d(zCount,dvIdx)      = 1;
            Dv2_d(zCount,dvIdx)      = -1;
            Dy1_d(zCount,zIdx)       = -maxL(r,1);
            Dy2_d(zCount,zIdx)       = maxL(r,1);
            b1_d(zCount)             = maxL(r,1);

            c_0(zCount,1)            = val;
            zCount                   = zCount+1;
        end
    end
end


% Consistency check
if (zCount-1)~=numAV
    error('Generation of auxiliary variables z failed!')
end


Dv_d    = [sparse(numAV,numDV), speye(numAV);
            Dv1_d, -speye(numAV);
            Dv2_d, speye(numAV)];
Dy_d    = [Dy1_d; Dy2_d; sparse(numAV,numIV)];
h       = [sparse(numAV,1); b1_d; sparse(numAV,1)];   


%% Assemble output matrices

d_Av = [-Av', sparse(numPV,numAV);
        -speye(numAV+numDV);
        Dv_d];
    
d_Ay = [sparse(numAV+numDV+numPV,numIV);
        Dy_d];

d_b = [-c;
       sparse(numAV+numDV,1);
       h];

d_c = [b;
       -c_0];
  

lb_d    = zeros(numAV+numIV+numDV,1);
ub_d    = [maxL;b1_d; ones(numIV,1)];
   

end