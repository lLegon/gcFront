function [designTable,currDesigns,currFitness]=newReduceDesigns(model, currDesigns, currFitness, FitnessFcn, tol, deleteGenes, delsToInds, maxreductionsize)
% After designs identified, remove genes/reactions that are not
% contributing to the performance of a design

alreadyTested=zeros(size(currDesigns,1),1);

% remove designs with same metrics but extra deletions
%[currDesigns,currFitness,alreadyTested]=removeSubsets(currDesigns,currFitness,alreadyTested);

% store all calculated designs
allFitness=currFitness;
allDesigns=currDesigns;

iteration=1;

% showStep=false;
% showTiming=false;

currDesigns=uint8(currDesigns);
% 
% [~,bb]=find(currDesigns);
% bb=unique(bb);
% 
% designChar=strrep(string(num2str(currDesigns(:,unique(bb)))),' ','');
% t3=table(currFitness,designChar);
% disp(t3)

while 1
    
    % make new designs by removing 1 deletion from existing designs
    reductionLevel=1;
    while 1
        
        [newDesigns,parentTracker,parentTrackerInt,parentsUsed]=removeOneDel(currDesigns,reductionLevel,alreadyTested);
        if ~isempty(newDesigns) || reductionLevel==4
            % stop trying to make new designs if new designs were
            % successfully made, or if reduction level cannot be raised
            % further
            break
        else
            % if no new designs could be made, raise the reduction level
            % and try making new designs again
            reductionLevel=reductionLevel+1;
        end
    end

    % end reduction if designs cannot be reduced further
    if isempty(newDesigns)
        fprintf('\nTerminating removal of redundant KOs- no further reduction possible')
        break
    else
        fprintf(['\nRemoval of redundant KOs round ',num2str(iteration),', Reduction strictness = ',num2str(reductionLevel),'/4, New designs: '])
    end

    % get fitness of new designs, and save these into allDesigns
    [newFitness,allDesigns,allFitness]=getDesignFitness(newDesigns,FitnessFcn,allDesigns,allFitness,tol);

    
    % find pareto front from these designs
    paretoFitness=getParetoDesigns([newFitness;currFitness]);
    

    % classify deletions based on how their loss affects fitness
    [newDesigns,newOnFront]=classifyDeletions(newDesigns,newFitness,paretoFitness,parentTracker,parentTrackerInt,currFitness(parentsUsed,:),tol);

    % make minimal designs with only deletions that are needed to stay on
    % pareto, and deletions that have same effect on fitness
    minimalDesigns=makeMinimalDesigns(newDesigns,newFitness,newOnFront,parentTracker,parentTrackerInt);

    
    if ~isempty(minimalDesigns)
        % get minimal design fitness
        fprintf(', New minimal designs: ')
        [minimalFitness,allDesigns,allFitness]=getDesignFitness(minimalDesigns,FitnessFcn,allDesigns,allFitness,tol);
        
        % Remove minimal set from other deletions with same parents
        
        % recalculate pareto with minimal designs
        paretoFitness=getParetoDesigns([paretoFitness;minimalFitness]);
        
%        designChar=strrep(string(num2str([newDesigns(:,unique(bb));minimalDesigns(:,unique(bb))])),' ','');
%        t=table([newOnFront;nan(size(minimalDesigns,1),1)],double(ismember([newFitness;minimalFitness],paretoFitness,'rows')),[newFitness;minimalFitness],designChar);
%         disp(t)
%         figure(f1)
%                     imagesc([newDesigns(:,bb);minimalDesigns(:,bb)])
%                     drawnow
        

        % store minimal designs in newDesigns
        newDesigns=[newDesigns;minimalDesigns];
        newFitness=[newFitness;minimalFitness];

        % find pareto front of designs after new designs created
        newOnFront=ismember(newFitness,paretoFitness,'rows');
    else
        fprintf(', New minimal designs: 0')
    end
    
    % find index of current designs that are on the pareto front
    currOnFront=ismember(currFitness,paretoFitness,'rows');
    
   % increase reduction level of any design that was reduced and is still
   % on the pareto front
    alreadyTested=alreadyTested(currOnFront);
    alreadyTested(alreadyTested<reductionLevel)=reductionLevel;
    
    % if there are no further deletions that will be made at a higher
    % reduction level, then increase reduction level further
    tempMat=currDesigns(currOnFront,:);
    tempMat(tempMat<=reductionLevel)=6;
    lowestNotReduced=min(tempMat,[],2)-1;
    alreadyTested=max([alreadyTested,lowestNotReduced],[],2);
    
    % construct a new currDesigns from the designs on the front
    currDesigns=[currDesigns(currOnFront,:);newDesigns(newOnFront,:)];
    currFitness=[currFitness(currOnFront,:);newFitness(newOnFront,:)];
    alreadyTested=[alreadyTested;zeros(sum(newOnFront),1)];
    
    % prevent further analysis of designs larger than maxreductionsize (doing this now in case
    % creation of minimal designs made something small enough for further
    % analysis
    if iteration==1
        alreadyTested(sum(currDesigns~=0,2)>maxreductionsize)=5;
    end

    % Remove designs with same metrics but extra deletions- don't waste time on
    % reducing them further
    [currDesigns,currFitness,alreadyTested]=removeSubsets(currDesigns,currFitness,alreadyTested);

%     designChar=strrep(string(num2str(currDesigns(:,unique(bb)))),' ','');
%     t2=table(currFitness,designChar,alreadyTested);
%     disp(t2)
%     figure(f2)
%     imagesc(currDesigns(:,bb))
% %     drawnow
%     disp(size(paretoFitness))
    
    iteration=iteration+1;
end

% create a new designTable that contains the newly calculated pareto
% front
currDesigns=currDesigns~=0;

[currDesigns,index]=unique(currDesigns,'rows');
currFitness=currFitness(index,:);

if ~deleteGenes
    designTable=table(strings(size(currDesigns,1),1),zeros(size(currDesigns,1),1),zeros(size(currDesigns,1),1),zeros(size(currDesigns,1),1),zeros(size(currDesigns,1),1),nan(size(currDesigns,1),1),'VariableNames',{'ReactionDeletions','NoOfDels','GrowthRate','ProductFlux','CouplingStrength','DistFromIdeal'});
else
    designTable=table(strings(size(currDesigns,1),1),zeros(size(currDesigns,1),1),zeros(size(currDesigns,1),1),zeros(size(currDesigns,1),1),zeros(size(currDesigns,1),1),nan(size(currDesigns,1),1),'VariableNames',{'GeneDeletions','NoOfDels','GrowthRate','ProductFlux','CouplingStrength','DistFromIdeal'});
end
for a=1:size(designTable,1)
    if ~deleteGenes
        designTable{a,1}=join(model.rxns(delsToInds(currDesigns(a,:))),' ');
    else
        designTable{a,1}=join(model.genes(delsToInds(currDesigns(a,:))),' ');
    end
    if deleteGenes
        designTable{a,2}=sum(currDesigns(a,:))+count(designTable{a,1},'+');
    else
        designTable{a,2}=sum(currDesigns(a,:));
    end
    
    designTable{a,3:5}=currFitness(a,:);
end

end

function [minimalDesigns]=makeMinimalDesigns(newDesigns,newFitness,newOnFront,parentTracker,parentTrackerInt)

% 1 minimal design which is all deletions that are not necessary to be on
% front
% also minimal designs consisting of all single deletions that give same
% fitness
minimalDesigns=uint8(zeros(size(newDesigns)));
i=1;

for a=1:length(parentTrackerInt)-1
    childDesigns=parentTrackerInt(a)+1:parentTrackerInt(a+1);
    frontChild=childDesigns(newOnFront(childDesigns));
    
    if isempty(frontChild) || length(frontChild)==1
        continue
    end
    
    % delete everything that could be lost without removing parent from
    % pareto
    minimalDesigns(i,:)=newDesigns(frontChild(1),:);
    minimalDesigns(i,parentTracker(frontChild))=0;
    i=i+1;
    
    % find single deletions where same point on pareto is achieved if they
    % are lost- try deleting them all at same time
    uniqueFitness=unique(newFitness(frontChild,:),'rows');
    for b=1:size(uniqueFitness,1)
        
        matchingFitness=childDesigns(ismember(newFitness(childDesigns,:),uniqueFitness(b,:),'rows'));
        
        if length(matchingFitness)==1
            % dont get this fitness by removing other deletions, so they
            % are distinct- no minimal design to make
            continue
        else
           minimalDesigns(i,:)=newDesigns(matchingFitness(1),:);
           minimalDesigns(i,parentTracker(matchingFitness))=0;
           i=i+1;
        end
        
    end
    
end

if i>1
    minimalDesigns=minimalDesigns(1:i-1,:);
else
    minimalDesigns=[];
end

end

function [newDesigns,onFront]=classifyDeletions(newDesigns,newFitness,paretoFitness,parentTracker,parentTrackerInt,parentFitness,tol)

% [~,bb]=find(newDesigns~=0);
% bb=unique(bb);
% subplot(1,2,1)
% imagesc(newDesigns(:,bb));

onFront=ismember(newFitness,paretoFitness,'rows');
uncoupled=newFitness(:,2)<tol;
notOnFront=and(~onFront,~uncoupled);

changeMatrix=false(size(newDesigns,1),4);

for a=1:length(parentTrackerInt)-1
    
    % skip if parent produced nothing on Pareto
    if sum(onFront(parentTrackerInt(a)+1:parentTrackerInt(a+1)))==0
        continue
    end
    
    % make matrix to store which child designs are related, and how each of
    % their deletions affected performance
    changeMatrix(:,:)=false;
    changeMatrix(parentTrackerInt(a)+1:parentTrackerInt(a+1),:)=true;
    
    % column 1: is child on front?
    % column 2: does removed deletion in this child affect parental design's fitness?
    % column 3: does removed deletion in this child allow coupling, but move design off pareto front?
    % column 4: does removed deletion in this child prevent coupling?
    changeMatrix(:,[1,3,4])=and(changeMatrix(:,[1,3,4]),[onFront,notOnFront,uncoupled]);
    if ismember(parentFitness(a,:),paretoFitness,'rows')
        changeMatrix(parentTrackerInt(a)+1:parentTrackerInt(a+1),2)=~ismember(newFitness(parentTrackerInt(a)+1:parentTrackerInt(a+1),:),parentFitness(a,:),'rows');
    end
    
    newDesigns(changeMatrix(:,1),parentTracker(changeMatrix(:,1)))=1;
    newDesigns(changeMatrix(:,1),parentTracker(changeMatrix(:,2)))=2;
    newDesigns(changeMatrix(:,1),parentTracker(changeMatrix(:,3)))=3;
    newDesigns(changeMatrix(:,1),parentTracker(changeMatrix(:,4)))=4;

end

% setting KOs to 1 or 2 can overwrite a KO being removed- so remake all the KOs
newDesigns( sub2ind( size(newDesigns),1:size(newDesigns,1),parentTracker') )=0;

% subplot(1,2,2)
% imagesc(newDesigns(:,unique(bb)))
% drawnow

end

function [paretoFitness]=getParetoDesigns(fitness)

% get unique fitness values so there will never be draws
fitness=unique(fitness,'rows');

% sort such that designs cannot be dominated by things below them
fitness=sortrows(fitness,[-1,-2,-3]);
dominated=false(size(fitness,1),1);

% loop over fitness and remove anything that is dominated by the existing
% design
for a=1:size(fitness,1)-1
    % if a design is not on Pareto front, then anything it dominates is also not
    % on front
    if dominated(a)
        continue
    end
    
    % designs that are not better than the current design in at least one
    % metric are dominated by it
    dominated(a+1:end)=or(dominated(a+1:end),sum(fitness(a+1:end,:)>fitness(a,:),2)==0);
    
end

% Pareto front = all non-dominated fitness values
paretoFitness=fitness(~dominated,:);

end

function [newFitness,allDesigns,allFitness]=getDesignFitness(newDesigns,FitnessFcn,allDesigns,allFitness,tol)
% check if design fitness has already been calculated
[dontTest,matchingRow]=ismember(newDesigns~=0,allDesigns,'rows');

% display the number of new designs to test
fprintf(num2str(sum(~dontTest)))


% retrieve stored fitness for previously tested designs
newFitness=zeros(size(newDesigns,1),3);
newFitness(dontTest,:)=allFitness(matchingRow(dontTest),:);

% calculate fitness for untested designs
newFitness(~dontTest,:)=-FitnessFcn(newDesigns(~dontTest,:)~=0);

% set fitness for uncoupled designs to 0 so they can't get onto pareto by
% having high growth
uncoupled=newFitness(:,2)<tol;
newFitness(uncoupled,:)=0;

[newUniqueDesigns,uniqueInd]=unique(newDesigns(~dontTest,:)~=0,'rows');
calculatedFitness=newFitness(~dontTest,:);

allDesigns=[allDesigns;newUniqueDesigns];
allFitness=[allFitness;calculatedFitness(uniqueInd,:)];

end

function [newDesigns,parentTracker,parentTrackerInt,parentsUsed]=removeOneDel(currDesigns,reductionLevel,alreadyTested)
% create matrix which contains creates a matrix of all designs that can be
% created by removing a single deletion from existing designs.
% Third parameter lets you avoid removing some deletions- 

% find number of new designs that can be created from existing designs
removableDels=currDesigns==reductionLevel;
removableDels(alreadyTested>=reductionLevel,:)=false;
nNewDesigns=sum(removableDels,2);
parentsUsed=nNewDesigns>0;

% create matrix to store new designs, and an array specifying their parent
newDesigns=uint8(zeros(sum(nNewDesigns),size(currDesigns,2)));
parentTracker=zeros(sum(nNewDesigns),1);

i=0;
for a=1:size(currDesigns,1)
    % fill matrix with all the designs that can be made from parent design
    newDesigns(i+1:i+nNewDesigns(a),:)=repmat(currDesigns(a,:),[nNewDesigns(a),1]);
    delsToRemove=find(removableDels(a,:)==1);
    for b=1:nNewDesigns(a)
        newDesigns(i+b,delsToRemove(b))=false;
    end
    % specify what was deleted
    parentTracker(i+1:i+nNewDesigns(a),1)=delsToRemove;
    % update index
    i=i+nNewDesigns(a);
end

% keep track of which deletions in parentTracker represent a new parent
parentTrackerInt=[0;cumsum(nNewDesigns)];
parentTrackerInt=unique(parentTrackerInt);

end

function [currDesigns,currFitnesses,alreadyTested]=removeSubsets(currDesigns,currFitnesses,alreadyTested)
    % find designs which have the same fitness but extra deletions when
    % compared to other designs
    
    % remove non-unique designs (note currDesigns is always above newDesigns, so newDesigns will be removed if they are identical to previous) 
    [currDesigns,uniqueIndex]=unique(currDesigns,'rows');
    currFitnesses=currFitnesses(uniqueIndex,:);
    alreadyTested=alreadyTested(uniqueIndex);
    % could have problems with removing a higher value of alreadyTested?
    
    % sort designs by metrics and number of deletions
    nDels=sum(currDesigns~=0,2);
    expandedMetrics=[currFitnesses,nDels];
    
    [expandedMetrics,index]=sortrows(expandedMetrics,[-1,-2,-3,4]);
    currDesigns=currDesigns(index,:);
    alreadyTested=alreadyTested(index);
    currFitnesses=expandedMetrics(:,1:3);
    nDels=expandedMetrics(:,4);
    
    
    toRemove=false(size(currDesigns,1),1);
    setEnd=0;
    for a=1:size(currDesigns,1)

        % skip if design has already been removed
        if toRemove(a)
            continue
        end
        
        % find where the start and end of designs with same fitness is
        if a>setEnd
            setStart=a;
            setEnd=a;
            while 1
                if setEnd==size(currDesigns,1) || ~isequal(currFitnesses(setStart,:),currFitnesses(setEnd+1,:))
                    break
                else
                    setEnd=setEnd+1;
                end
            end
        end
        
        % skip if no other designs have same fitness
        if setStart==setEnd
            continue
        end
        
        % create index of designs with identical metrics
        matchIndVals=setStart:setEnd;
        
        % only consider designs with more deletions than this one
        moreDels=nDels(matchIndVals)>nDels(a);
        matchIndVals=matchIndVals(moreDels);
        if isempty(matchIndVals)
            continue
        end
        
        % remove designs with more deletions than this one
        toRemove(matchIndVals)=true;
%         
%         % find any of these designs with all the same deletions
%         matchingDels=tempDesigns(matchIndVals,tempDesigns(a,:));
%         toRemove(matchIndVals(sum(matchingDels>0,2)==size(matchingDels,2)))=true;
        
    end
    
    currDesigns=currDesigns(~toRemove,:);
    currFitnesses=currFitnesses(~toRemove,:);
    alreadyTested=alreadyTested(~toRemove);
    
end