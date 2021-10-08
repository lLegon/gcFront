function [designTable]=reduceDesigns(model, currDesigns, currFitnesses, FitnessFcn, tol, deleteGenes, delsToInds)
% After designs identified, remove genes/reactions that are not 
% contributing to the performance of a design

    completeDesigns=false(size(currDesigns,1),1);
    
    allTestedDesigns=currDesigns;
    
    while 1
        
        nNewDesigns=sum(currDesigns,2);
        nNewDesigns(completeDesigns)=0;
        if sum(nNewDesigns)==0
            % end loop if every design has already been reduced
            break
        end
        
        % create matrix to store the designs to be tested
        newDesigns=false(sum(nNewDesigns),size(currDesigns,2));
        
        
        % loop through existing designs and create new designs by removing
        % one deletion from existing
        dInd=1;
        
        for a=1:size(currDesigns,1)
            
            
            if completeDesigns(a)==true
                % skip if design has already been reduced
                continue
            end
            
            % store the current design in the newDesigns matrix
            newDesigns(dInd:dInd+nNewDesigns(a)-1,:)=repmat(currDesigns(a,:),[nNewDesigns(a),1]);
            
            % remove one deletion from each of the designs in newDesigns
            deletions=find(currDesigns(a,:));
            for b=1:length(deletions)
                newDesigns(dInd+b-1,deletions(b))=false;
            end
            
            % update index
            dInd=dInd+nNewDesigns(a);
            
        end
        
        % remove any non-unique designs
        newDesigns=unique(newDesigns,'rows');
        
        % remove any designs that are already in allTestedDesigns
        
        alreadyTested=ismember(newDesigns,allTestedDesigns,'rows');
        
        newDesigns=newDesigns(~alreadyTested,:);
        
        % stop if there are no designs left that haven't been tested already
        if isempty(newDesigns)
            break
        end
        
        
        
        % find fitnesses using fitnessFcn
        
        newFitnesses=-feval(FitnessFcn,newDesigns);
        
        
        % keep any of the designs that are GC
        gcInds=newFitnesses(:,2)>tol;
        
        % combine old and new designs
        tempDesigns=[currDesigns;newDesigns(gcInds,:)];
        tempFitnesses=[currFitnesses;newFitnesses(gcInds,:)];
        tempCompleted=[true(size(currDesigns,1),1);false(sum(gcInds),1)];
        allTestedDesigns=[allTestedDesigns;newDesigns]; %#ok<AGROW>
        
        % sort designs by fitness so designs with the same metrics can be
        % assessed together
        [tempFitnesses,sortInd]=sortrows(tempFitnesses,[1,2,3]);
        tempDesigns=tempDesigns(sortInd,:);
        tempCompleted=tempCompleted(sortInd);
        
        
        
        % calculate pareto front
        onFront=false(size(tempDesigns,1),1);
        for a=1:size(tempFitnesses,1)
            % if previous design has same fitness scores, then this design
            % is on the front if the previous one was
            if a~=1 && isequal(tempFitnesses(a,:),tempFitnesses(a-1,:))
                onFront(a)=onFront(a-1);
                continue
            end
            
            % to find points that dominate current point, search for points
            % that are equal or better at everything, and then remove any
            % points with exactly the same score for everything
            dominant=and(sum(tempFitnesses>=tempFitnesses(a,:),2)==size(tempFitnesses,2),~ismember(tempFitnesses,tempFitnesses(a,:),'rows'));
            
            % if no dominant point exists, then this point is on the front
            if sum(dominant)==0
                onFront(a)=true;
            end
            
        end
        
        % find and remove designs if they have the same metrics and
        % deletions as a second design, but also have a few more deletions, since
        % further reduction will just lead to the second design,
        % so is pointless
        uFit=unique(tempFitnesses(onFront,:),'rows');
        minimalDesigns=false(size(uFit,1),size(tempDesigns,2));
        for b=1:size(uFit,1)
            % find all designs that give a particular metric
            matchingInds=find(ismember(tempFitnesses,uFit(b,:),'rows'));
            if length(matchingInds)==1
                % if only one design gives this metric, then continue
                continue
            end
            % loop over designs, and find if there is anything that has all
            % the same deletions, but also some extras that don't do
            % anything
            
            redundantDels=false(length(matchingInds),1);
            for c=1:length(matchingInds)
                if redundantDels(c)==true
                    % no point in analysing a design already shown to be
                    % redundant, as any matching design with even more
                    % redundancy will have already been identified when
                    % this design was identified
                    continue
                end
                
                % find designs with same deletions as the current design,
                % and label them as containing redundant dels
                matchingDels=sum(tempDesigns(matchingInds,tempDesigns(matchingInds(c),:)),2);
                redundantDels=or(redundantDels,matchingDels==matchingDels(c));
                redundantDels(c)=false;
            end
            % remove redundant KO combinations from Pareto
            if sum(redundantDels)~=0
                onFront(matchingInds(redundantDels))=false;
            end
            % if there are multiple non-redundant designs, create a minimal
            % design that only contains deletions common to all of
            % them
            if sum(~redundantDels)>1
                minimalDesigns(b,:)=sum(tempDesigns(matchingInds(~redundantDels),:),1)==sum(~redundantDels);
            end
        end
        
        % remove non-unique minimal designs
        minimalDesigns=unique(minimalDesigns,'rows');
        % remove any designs that are already in allTestedDesigns
        alreadyTested2=ismember(minimalDesigns,allTestedDesigns,'rows');
        minimalDesigns=minimalDesigns(~alreadyTested2,:);
        
        % stop if there are no designs left that haven't been tested already
        
        if ~isempty(minimalDesigns)
            
            minFitnesses=-feval(FitnessFcn,minimalDesigns);
            gcInds2=minFitnesses(:,2)>tol;
            
            % add minimal designs to tempDesigns
            tempDesigns=[tempDesigns(onFront,:);minimalDesigns(gcInds2,:)];
            tempFitnesses=[tempFitnesses(onFront,:);minFitnesses(gcInds2,:)];
            tempCompleted=[tempCompleted(onFront,:);false(sum(gcInds),1)];
            allTestedDesigns=[allTestedDesigns;minimalDesigns]; %#ok<AGROW>
            
            % recalculate pareto and remove redundant combinations
            
            
            % sort designs by fitness so designs with the same metrics can be
            % assessed together
            [tempFitnesses,sortInd]=sortrows(tempFitnesses,[1,2,3]);
            tempDesigns=tempDesigns(sortInd,:);
            tempCompleted=tempCompleted(sortInd);
            
            
            
            % calculate pareto front
            onFront=false(size(tempDesigns,1),1);
            for a=1:size(tempFitnesses,1)
                % if previous design has same fitness scores, then this design
                % is on the front if the previous one was
                if a~=1 && isequal(tempFitnesses(a,:),tempFitnesses(a-1,:))
                    onFront(a)=onFront(a-1);
                    continue
                end
                
                % to find points that dominate current point, search for points
                % that are equal or better at everything, and then remove any
                % points with exactly the same score for everything
                dominant=and(sum(tempFitnesses>=tempFitnesses(a,:),2)==size(tempFitnesses,2),~ismember(tempFitnesses,tempFitnesses(a,:),'rows'));
                
                % if no dominant point exists, then this point is on the front
                if sum(dominant)==0
                    onFront(a)=true;
                end
                
            end
            
            % find and remove designs if they have the same metrics and
            % deletions as a second design, but also have a few more deletions, since
            % further reduction will just lead to the second design,
            % so is pointless
            uFit=unique(tempFitnesses(onFront,:),'rows');
            for b=1:size(uFit,1)
                % find all designs that give a particular metric
                matchingInds=find(ismember(tempFitnesses,uFit(b,:),'rows'));
                if length(matchingInds)==1
                    % if only one design gives this metric, then continue
                    continue
                end
                % loop over designs, and find if there is anything that has all
                % the same deletions, but also some extras that don't do
                % anything
                
                redundantDels=false(length(matchingInds),1);
                for c=1:length(matchingInds)
                    if redundantDels(c)==true
                        % no point in analysing a design already shown to be
                        % redundant, as any matching design with even more
                        % redundancy will have already been identified when
                        % this design was identified
                        continue
                    end
                    
                    % find designs with same deletions as the current design,
                    % and label them as containing redundant dels
                    matchingDels=sum(tempDesigns(matchingInds,tempDesigns(matchingInds(c),:)),2);
                    redundantDels=or(redundantDels,matchingDels==matchingDels(c));
                    redundantDels(c)=false;
                end
                % remove redundant KO combinations from Pareto
                if sum(redundantDels)~=0
                    onFront(matchingInds(redundantDels))=false;
                end
                
            end
            
        end
        
        currDesigns=tempDesigns(onFront,:);
        currFitnesses=tempFitnesses(onFront,:);
        completeDesigns=tempCompleted(onFront,:);
        
       
        
    end
    
    % create a new designTable that contains the newly calculated pareto
    % front
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
        
        designTable{a,3:5}=currFitnesses(a,:);
    end

end