function [Av,Ay,b,numKORxns,nonKnock,nonKnock_s] = genMIPP(nRxns,nRxns_s,...
                    S_s,ub_s,lb_s,nMets,notknockables,mapRxnsRed,B)

% Generate matrices of general MIP problem
Av      = [];
Ay      = [];
b       = [];
I       = eye(nRxns_s);

% Match non-knockables and optional external fluxes to split network
[~,nonKnock]        = find(mapRxnsRed(notknockables,:));
[~,nonKnock_s]      = find(B(nonKnock,:));


numNonKO            = length(nonKnock);

numKORxns   = nRxns-numNonKO;


% Special constraints (substrate uptake, ATP maintenance etc.) -> include i
% general constraints!

% Constraints non-knockable reactions
Av      = [Av;I(nonKnock_s,:);-I(nonKnock_s,:)];
Ay      = [Ay;zeros(2*size(nonKnock_s,1),numKORxns)];
b       = [b;ub_s(nonKnock_s);-lb_s(nonKnock_s)];

% Equality constraints S*v=0
Av      = [Av;S_s;-S_s];
Ay      = [Ay;zeros(nMets,numKORxns);zeros(nMets,numKORxns)];
b       = [b;zeros(nMets,1);zeros(nMets,1)];


% KO assignments + flux boundaries (implying positive and knockable fluxes)
Iknockable                  = I;
Iknockable(nonKnock_s,:)    = [];
ubKnockable                 = ub_s;
lbKnockable                 = lb_s;
ubKnockable(nonKnock_s)     = [];
lbKnockable(nonKnock_s)     = [];
Bknockable                  = B;
Bknockable(:,nonKnock_s)    = [];
Bknockable(nonKnock,:)      = [];

Ay1     = Bknockable*diag(ubKnockable);

Ay2     = Bknockable*diag(lbKnockable);

Av      = [Av;Iknockable;-Iknockable];          
Ay      = [Ay;-Ay1';Ay2'];
b       = [b; zeros(2*(nRxns_s-size(nonKnock_s,1)),1)];



end
