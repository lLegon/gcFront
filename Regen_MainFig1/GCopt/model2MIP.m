function [A_j, b_j, lb_j, ub_j, cO, vSize, zSize, lSize, ySize, B_flux, mapRxnsRed,B, nonKnock, model_s,results_wt] =...  
                model2MIP(model,notKORxns,rxnObjI,rxnObjO,subsRxn,maxKO,reductionFlag,solver,probOpts,genOpts)


fprintf('--> Transform Model to MILP \n')

% Inner objective (normally max biomass)
c                   = zeros(size(model.rxns,1),1);
c(strcmp(model.rxns,rxnObjI))       = 1;

% Outer objective (OptKnock Style)
cO          = zeros(size(model.rxns,1),1);
cO(strcmp(model.rxns,rxnObjO))      = probOpts.sense;

% Extract rxn numbers of not-knockable reactions
if isempty(notKORxns)
    notknockables   = [];
else
    [~,notknockables]   = ismember(notKORxns,model.rxns);
end


%% Load model data
S       = model.S;
lb      = model.lb;
ub      = model.ub;
rev     = model.rev;
rxns    = model.rxns;
mets    = model.mets;

nRxns   = size(model.rxns,1);
nMets   = size(model.mets,1);


%% Determine product formation rate at optimal productivity
if isempty(probOpts.fixMu)
    solSpaceEdge    = prodSpaceCheck(model,rxnObjI,subsRxn,rxnObjO,probOpts.sense);

    results_wt      = solSpaceEdge;

    [maxProd,I]     = max(solSpaceEdge.prodR);
    disp(['Maximal productivity: ',num2str(maxProd),' mol/h'])
    disp(['Optimal production rate: ',num2str(solSpaceEdge.yieldR(I)),' mol/g/h'])
    disp(['Optimal Growth-Rate: ',num2str(solSpaceEdge.muR(I)),' /h'])
    optObjFlux      = solSpaceEdge.yieldR(I);
    optBmFlux       = solSpaceEdge.muR(I);

    results_wt.optObjFlux   = optObjFlux;
    results_wt.optBmFlux    = optBmFlux;
    results_wt.maxProd      = maxProd;
else
    solSpaceEdge    = prodSpaceCheck(model,rxnObjI,subsRxn,rxnObjO,probOpts.sense);
    results_wt      = solSpaceEdge;
    
    optBmFlux               = probOpts.fixMu;
    
    % get optimal production rate
    model_bm    = changeRxnBounds(model,rxnObjI,optBmFlux,'b');
    model_bm    = changeObjective(model_bm,rxnObjO);
    if probOpts.sense == 1
        FBASol      = optimizeCbModel(model_bm,'max');
    elseif probOpts.sense == -1
        FBASol      = optimizeCbModel(model_bm,'min');
    end
    
    results_wt.optObjFlux   = FBASol.f;
    results_wt.optBmFlux    = optBmFlux;
end


%% Compress network and target space
fprintf('\t ...Compress metabolic network \n')

% Reduce target space
notknockables =     redTargetSpace(model, notknockables, genOpts.modelType);


switch reductionFlag
    case 0  % No reduction
        mapRxnsRed      = eye(nRxns);

        % Transpose notknockable vector
        notknockables   = notknockables';
        
        
        model_s.reacidx         = eye(nRxns);
        model_s.metidx          = eye(nMets); 
        model_s.model           = model;
        
    case 1  % Only remove blocked reactions
        
        
        % FVA to identify blocked reactions
        minFlux     = zeros(length(model.rxns),1);
        maxFlux     = minFlux;
        parfor i=1:length(model.rxns)
            globalLPSolverVariable('gurobi');
            par_model   = model;
            % maximization
            par_model   = changeObjective(par_model,par_model.rxns(i));
            sol     = optimizeCbModel(par_model,'max');
            maxFlux(i)  = sol.x(i);
            % minimization
            sol     = optimizeCbModel(par_model,'min');
            minFlux(i)  = sol.x(i);
        end
            
      
        blockedRxns         = find(((minFlux==0) & (maxFlux==0))); %| ((lb==0) & (ub==0)));
        numBlockedRxns      = size(blockedRxns,1);
        reacidx             = eye(nRxns);   % Maps reactions of reduced network to original network
        metidx              = eye(nMets);   % Maps metabolites of reduced network to original network

        % Check removed reactions with objective function(s)
        objRxns         = find(c);
        for i=1:length(objRxns)
            if ~isempty(find(blockedRxns==objRxns(i),1))
                error('Considered reaction in objective function is blocked!')
            end
        end


        % remove reaction by constraining fluxes to zero
        model_s.model   = model;
        model_s.model   = changeRxnBounds(model_s.model,model_s.model.rxns(blockedRxns),0,'b');
        
        
        nRxns       = size(rxns,1);
        mapRxnsRed  = eye(nRxns);
        
        % Display blocked reactions
        fprintf(['\t\t Blocked Reactions: ',num2str(numBlockedRxns),'\n'])
%         disp(model.rxnNames(blockedRxns))

        % Remove unused metabolites
        delMet                  = any(S,2)==0;
        S(delMet,:)             = [];
        mets(delMet)            = [];
        metidx(:,delMet)        = [];
        
        nMets   = size(mets,1);
        
        % Display removed metabolites
        fprintf(['\t\t Removed metabolites: ',num2str(sum(delMet)),'\n'])
%         disp(model.metNames(delMet))
        
           
        % (Re-)Assign non-knockables
        notKnockCompr   = any(reacidx(notknockables,:),1);
        notknockables   = find(notKnockCompr);
        
 
        model_s.reacidx         = reacidx;
        model_s.metidx          = metidx; 
        
             
end


% determine new network size
revRxns     = find(rev);
irrevRxns   = find(rev==0);
nRev        = size(revRxns,1);
nIrrev      = size(irrevRxns,1);
nRxns_s     = nRev+nRxns;


%% Split reversible reactions
S_s         = [S, zeros(nMets,nRev)];    % stoichiometric matrix after splitting
lb_s        = [lb; zeros(nRev,1)];
ub_s        = [ub; zeros(nRev,1)];
c           = [c; zeros(nRev,1)];
cO          = [cO; zeros(nRev,1)];
B           = [eye(nRxns),zeros(nRxns,nRev)];       % Mapping matrix B
B_flux      = B;
objRxns_r   = find(c);
objRxnsO_r  = find(cO);


for i=1:nRev
    S_s(:,nRxns+i)   = -S(:,revRxns(i));
    ub_s(nRxns+i)    = -lb(revRxns(i));
    lb_s(revRxns(i)) = 0;
    B(revRxns(i),nRxns+i)       = 1;
    B_flux(revRxns(i),nRxns+i)  = -1;
    
    if any(objRxns_r==revRxns(i))
        c(nRxns+i,1)    = -1;           % If objective reaction is reversible, its negative reverse flux must be maximised
    end
    if any(objRxnsO_r==revRxns(i))
        cO(nRxns+i,1)    = -1;           % If objective reaction is reversible, its negative reverse flux must be maximised
    end    
     
end


%% Determine essential reactions and exclude from target space
[eRxns,~]       = essentialRxns(S_s,lb_s,ub_s,c,B);
fprintf(['\t\t Essential reactions: ',num2str(sum(eRxns)),'\n'])
notknockables   = unique([notknockables,find(eRxns)']);

fprintf(['\t\t',num2str(size(S,2)-length(notknockables)),' target reactions out of ',num2str(length(model.rxns)),'\n']) 


%% Fix growth rate to optimal (productivity) value
ub_s(objRxns_r)         = optBmFlux;
lb_s(objRxns_r)         = optBmFlux;



%% Generate matrices of a general MIP problem
[Av,Ay,b,numKORxns,nonKnock,nonKnock_s]   = genMIPP(nRxns,nRxns_s,...
                    S_s,ub_s,lb_s,nMets,notknockables,mapRxnsRed,B);


% All entries in Ay1/Ay2
[indLR, indLC]  = find(Ay);
zSize   = size(indLR,1);      % Number of auxiliary variables z


lSize   = size(Av,1);                           % Number of dual variables lambda
ySize   = numKORxns;                            % Number of integer variables (knockable reactions before split)
vSize   = nRxns_s;                              % Number of reactions after the split


% Define production rate as objective (to be minimised)
cOpt        = -cO;

%% Dualise
[d_Av, d_Ay, d_b, d_c, lb_d, ub_d] = dualise(Av, Ay, b, cOpt, solver, indLR, indLC, genOpts);


%% Join dual and primal problem

d_vSize         = size(d_Av,2);
numPrimConstr   = size(Av,1);
numDualConstr   = size(d_Av,1);


% Assemble matrices of the final optimisation problem
A_j     = [cOpt', -d_c', sparse(1,ySize);                         % equality between dual and primal objective
            Av, sparse(numPrimConstr,d_vSize), Ay;                      % Stoichiometric constraints / flux boundaries of primal problem
            sparse(numDualConstr,vSize), d_Av, d_Ay;                    % Constraints of the dual problem
            sparse(1,vSize+d_vSize),-ones(1,ySize)];   % Max KO number constraint
b_j     = [0;
            b;
            d_b;
            maxKO-ySize];

lb_j    = [lb_s; lb_d];

ub_j    = [ub_s; ub_d];
    
    

%% Save reduced and split model network
model_s.S_s         = S_s;
model_s.ub_s        = ub_s;
model_s.lb_s        = lb_s;
model_s.c_s         = c;
model_s.cO_s        = cO;
model_s.nonKnock_s  = nonKnock_s;
model_s.nonKnock    = nonKnock;
model_s.bmRxn       = rxnObjI;
model_s.targetRxn   = rxnObjO;
model_s.subsRxn     = subsRxn;
model_s.B           = B;

end


function globalLPSolverVariable(solver)
% workaround to allow for the use of COBRA toolbox function
% "optimizeCbModel" within parfor loop
global CBT_LP_SOLVER
CBT_LP_SOLVER   = solver;
end

