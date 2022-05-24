function [fitness2] = myFitness3obj(delInds, poolModel, delsToInds, targetInd, growthInd, nRxns, tol, minGrowth, delta, defaultFitness, tiltVal, deleteGenes, geneRules, geneMat)
% calculate fitness of a list of deletions
% delInds: matrix of deletions to test
% poolModel: metabolic model stored as a parallel.pool.Constant object. If the parallel pool has closed, this function will fail
% delsToInds: vector that links reaction/gene deletions stored in x to the index of this reaction/gene in the model
% targetInd: index of the reaction that is being optimised
% growthInd: index of the biomass reaction
% nRxns: total number of reactions in the model
% tol: rounding tolerance used to control mathematical errors
% minGrowth: minimum growth threshold
% delta: how much product synthesis should be increased by when checking shadow price for uncoupled designs
% defaultFitness: fitness value to be assigned if deletion makes model infeasible
% tiltVal: the value used to tilt the objective vector
% deleteGenes: whether genes or reactions are being deleted
% geneRules: the model's reaction-gene rules
% geneMat: the model's reaction-gene matrix
% note that the fitness function will assign uncoupled designs a growth of
% zero if coupled designs have been found, to ensure that uncoupled designs
% cannot be on the pareto



% removing duplicates from x
[delInds,~,shortToLong]=unique(delInds,'rows');


if ~deleteGenes
    % convert deletions from position in selected reaction list -> position in model
    delInds2=false(size(delInds,1),nRxns);
    delInds2(:,delsToInds)=delInds;
else
    % convert gene deletions from position in selected gene list -> position in model
    delIndsTemp=false(size(delInds,1), size(geneMat,2));
    delIndsTemp(:,delsToInds)=delInds;
    
    %tic
    % convert gene deletions into rxn deletions
    delIndsTemp2=false(size(delInds,1),nRxns);
    for a=1:size(delIndsTemp2,1)
        x=~delIndsTemp(a,:);
        [assocRxns,~]=find(geneMat(:,delIndsTemp(a,:)));
        assocRxns=unique(assocRxns);
        for b=1:length(assocRxns)
            
            delIndsTemp2(a,assocRxns(b))=~eval(geneRules{assocRxns(b)});
        end
    end

    
    [delInds2,~,shortToLong2]=unique(delIndsTemp2,'rows');
    
    %toc
end


% creating fitness vector
nCombinations=size(delInds2,1);
fitness=zeros(nCombinations,3);
fitness(:,3)=defaultFitness;


environment=getEnvironment();

%tic
parfor a=1:nCombinations
   
    restoreEnvironment(environment,0);
    
    % generating KO model
    
    KOmodel=poolModel.Value;
%      KOmodel=poolModel;
     
    KOmodel.lb(delInds2(a,:))=0;
    KOmodel.ub(delInds2(a,:))=0;
    
    s1=optimizeCbModel2(KOmodel);
    
    if s1.stat~=1 || s1.x(growthInd)<=minGrowth+tol
        continue
    end
    
    if s1.x(targetInd)<tol
        % not growth coupled- applying first part of fitness function
        
        if s1.w(targetInd)==0%<=tiltVal+tol
            % reduced cost can sometimes become 0 even if it shouldn't be?
            % Sometimes happens if product synthesis is impossible, unclear
            % why it happens for other cases
            % if value is zero, should double check shadow price
            
            KOmodel.lb(targetInd)=delta;
            s2=optimizeCbModel2(KOmodel);
            
            % if increasing product synthesis makes model infeasible, then
            % coupling strength should be set to default
            if s2.stat~=1
                fitness(a,:)= [round(s1.x(growthInd)/tol)*tol,0,defaultFitness];
                continue
            end
            
            couplingStrength=(round((s1.x(growthInd)-s2.x(growthInd))/tol)*tol)/-delta;
            
            % if no product synthesis, coupling strength cannot be greater
            % than 0. This section is just in case mathematical errors
            % result in the coupling strength calculation giving a number
            % greater than 0.
            if couplingStrength>0
                couplingStrength=0;
            end

            if s2.stat==1 && couplingStrength > defaultFitness
                fitness(a,:)= [round(s1.x(growthInd)/tol)*tol,0,couplingStrength];

            end
        else
            
            if -round(s1.w(targetInd)/tol)*tol<0
                fitness(a,:)=[round(s1.x(growthInd)/tol)*tol,0,-round(s1.w(targetInd)/tol)*tol];
            else
                fitness(a,:)=[round(s1.x(growthInd)/tol)*tol,0,0];
            end
            
        end
        
        
    else
        growth=round(s1.x(growthInd)/tol)*tol;
        product=round(s1.x(targetInd)/tol)*tol;

        % minimising product synthesis while maximising growth
        KOmodel.lb(targetInd)=0;
        KOmodel.c(targetInd)=-1;
        KOmodel.c(growthInd)=tiltVal;
        
        s2=optimizeCbModel2(KOmodel);
        
        if s2.x(targetInd)>tol
            % ie product synthesis is necessary even if growth doesn't
            % occur
            couplingstrength=min( (round(s2.x(targetInd)/tol)*tol)/product, 1) + 1;
        else
            % product synthesis can be 0
            couplingstrength=(growth - round(s2.x(growthInd)/tol)*tol)/growth;
        end
        

        fitness(a,:)=[growth,product,couplingstrength];
        
        

            
     end
    
    
end
%toc

% if uncoupled designs are assigned a growth fitness score, then some of
% them will exist on the pareto. Since this is a waste, uncoupled designs
% will be assigned a growth of 0 once a coupled design has been discovered

if max(fitness(:,2))>tol
   fitness(fitness(:,2)<tol,1)=0;
end


% GA is trying to find a minimum, while we want the maximum, so use negative of calculated fitness
fitness=-fitness;

% converting fitness metric so it corresponds to original matrix
if ~deleteGenes
    fitness2=fitness(shortToLong,:);
else
    fitnessTemp=fitness(shortToLong2,:);
    fitness2=fitnessTemp(shortToLong,:);
end



end
