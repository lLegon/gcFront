% This script will find reaction knockouts that couple succinate to growth in the E. coli core model

% load model
model=readCbModel('e_coli_core.mat'); % either a structure containing a COBRA-toolbox compatible metabolic model, or the address of a file containing a COBRA-toolbox compatible metabolic model
% set glucose supply to 10 mmol/gDW/h
model=changeRxnBounds(model,"EX_glc__D_e",-10,'l'); % changing lower bound of EC_glc__D_e to -10
save('e_coli_core.mat','model'); % saving change in bounds

% call gcFront (options will now have to be entered manually into text boxes)
[designTable,algParams,reducedModel]=gcFront;




% % To supply gcFront's options from within a script instead of entering manually,
% % place a % sign at the start of line 10, and then use the "uncomment"
% % option on all the code after this line.
% 
% % set succinate as the metabolite to be coupled
% target='succinate'; % the reaction/metabolite in the model to be growth-coupled. if a metabolite name is supplied, then the exchange reaction for this metabolite will be selected
% 
% % set model options so gurobi is used to solve LP problems, and so knockout
% % combinations with more than 5 knockouts will not be considered
% options=struct;
% options.solver='gurobi';           % Character array of the name of the LP solver used for solving FBA problems. Default  =  currently set LP solver. If no solver set, will use whatever solver the function initCobraToolbox sets 
% options.maxknockouts=5;            % Double specifying the maximum number of knockouts that a design may contain. Default  =  inf (i.e. no limit on the number of deletions allowed)
% 
% % Other options that are not being set here:
% 
% options.biomassrxn='';           % Character array specifying the reaction that represents cell growth. Default  =  current objective reaction in model
% options.tol=10^-8;               % Double specifying tolerance to mathematical errors. Flux differences smaller than this will be ignored. Default  =  10^-8
% options.tiltval=10^-4;           % Double specifying the coefficient used for objective value tilting. Default  =  10^-4
% options.shiftval=10^-5;          % Double specifying the size of the change in growth rate used to calculate shadow price of uncoupled designs. Default  =  10^-5
% options.skipreduction=false;     % Logical that determines if algorithm should skip deletion of inactive reactions and pooling of unbranched pathways/redundant genes. Default  =  false (i.e. carry out reduction)
% options.mingrowth=10^-3;         % Double specifying the minimum growth threshold- do not consider any deletion that lower growth below this value. Default  =  10^-3
% options.minprod=10^-3;           % Double specifying the minimum product threshold- do not consider any deletion that lowers maxmimum product synthesis below this value. Default  =  10^-3
% options.removeredundancy=true;   % Logical that detemines if any redundant deletions should be removed from the designs that the GA identifies. Default  =  true
% options.newredundantremoval=true;% Logical specifying whether to use a new, faster method for removing redundant KOs, or to use the original version that was used to obtain the data for the gcFront paper.
% options.maxreductionsize=15;     % Double specifying the maximum size of design that should be fully explored for removal of redundant KOs
% options.saveresults=true;        % Logical that determines if results and algorithm parameters shoudl be saved. Default  =  true
% options.deletegenes=false;       % Logical that determines if algorithm should test gene knockouts or reaction knockouts. Default  =  false (i.e. algorithm will knock out reactions)
% options.ignorelistrxns={};       % cell array of reactions that user does not want to be knocked knocked out (e.g. reactions that are associated with experimentally inaccessible genes). Default  =  {} (i.e. no reactions)
% options.ignorelistgenes={};      % cell array of genes that should not be knocked out (e.g. essential genes). If reaction knockouts are being tested, reactions that can't be knocked out without knocking out these genes will be ignored. Default  =  {} (i.e. no genes).
% options.dontkoess=true;          % Logical to determine whether reactions that can only be knocked out if an essential gene is knocked out are ignored. Default  =  true (i.e do not consider these reactions)
% options.onlykogeneassoc=true;    % Logical to determine if reactions that are not gene associated should be ignored. Default  =  true (i.e. only consider deletion of gene associated reactions)
% options.mutationrate=1;          % Double controlling the mean number of mutations that the GA mutation function introduces. Default  = 1
% options.popsize=200;             % Double controlling the number of individuals in the GA population. Default  = 200
% options.genlimit=10000;          % Double specifying the maximum number of generations the GA will run for. Default  = 10000
% options.timelimit=60*60*24;      % Double specifying the maximum length of time (in seconds) that the GA will run for. Default  =  86400 (i.e. 1 day)
% options.fitnesslimit=inf;        % Double. The algorithm will terminate if a design with a product synthesis rate higher than this is discovered. Default = inf (i.e. this condition will not stop the algorithm)
% options.spreadchangelimit=10^-4; % Double specifying how low the change in spread must be before the algorithm terminates. For further information see MATLAB's 'gamultiobj Algorithm' documentation. Default  =  10^-4 
% options.stallgenlimit=[];        % Double specifying the number of generations that change in spread must be below spreadchangelimit before the algorithm terminates. For further information see MATLAB's 'gamultiobj Algorithm' documentation. Default  =  number of generations (i.e. the algorithm will never terminate due to spread being too low)
% options.plotinterval=1;          % Double specifying the number of generations that must pass before the algorithm updates the plot of the current designs. Default  =  1 (i.e. updates every generation)
% [designTable,algParams,reducedModel]=gcFront(model,target,options);
% 

