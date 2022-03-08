function pop=myCreation(delSetLength,~,options, mutationRate, maxKnockouts, geneCounter)
% creates the initial population for multiobjective GA

% create population according to mutation rate
popRands=rand(sum(options.PopulationSize-size(options.InitialPopulation,1)),delSetLength);
if size(options.InitialPopulation)==0
    pop=popRands<mutationRate;
else
    pop=[options.InitialPopulation;popRands<mutationRate];
end

if maxKnockouts~=inf
    % check population does not have too many knockouts
    if isempty(geneCounter)
        % rxn KOs- just sum the number of KOs in vector
        excessMutations=sum(pop,2)-maxKnockouts;
        for a=1:size(pop,1)
            if excessMutations(a)>0 % has more mutations than is allowed
                
                dels=find(pop(a,:));
                delsToRemove=dels(randperm(length(dels),excessMutations(a)));
                pop(a,delsToRemove)=false;
                
            end
        end
    else
        excessMutations=pop*geneCounter-maxKnockouts;
        for a=1:size(pop,1)
            if excessMutations(a)>0
                % find deletions
                dels=find(pop(a,:));
                delsToRemove=false(size(dels));
                delCounter=geneCounter(pop(a,:));
                % remove deletions until less than or equal to maxKnockouts
                while 1
                    delsToRemove(randi(length(delsToRemove),1))=true;
                    if sum(delCounter(~delsToRemove))<=maxKnockouts
                        break
                    end
                end
                % remove the specified deletions from the population
                pop(a,dels(delsToRemove))=false;
            end
        end
    end
end


end