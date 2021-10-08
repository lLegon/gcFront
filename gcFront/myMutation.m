function kids=myMutation(parents,options,GenomeLength,FitnessFcn,state,score,population,mutationRate,maxKnockouts,geneCounter)
% create random mutations in parents for multiobjective GA

kids=population(parents,:);
randInds=rand(size(kids))<mutationRate;
kids(randInds)=~kids(randInds);

% check if children have too many KOs

if maxKnockouts~=inf
    for a=1:size(kids,1)
        if isempty(geneCounter) 
            if sum(kids(a,:),2)>maxKnockouts
                % rxn KOs- remove knockouts so number of knockouts matches maximum number of
                % knockouts
                koInds=find(kids(a,:));
                removeInds=randperm(length(koInds),length(koInds)-maxKnockouts);
                kids(a,koInds(removeInds))=false;
            end
        elseif sum(geneCounter(kids(a,:)))>maxKnockouts
            % gene KOs- remove knockouts while taking into consideration
            % how many gene KOs each position entails
            dels=find(kids(a,:));
            delsToRemove=false(size(dels));
            delCounter=geneCounter(kids(a,:));
            % remove deletions until less than or equal to maxKnockouts
            while 1
                delsToRemove(randi(length(delsToRemove),1))=true;
                if sum(delCounter(~delsToRemove))<=maxKnockouts
                    break
                end
            end
            % remove the specified deletions from the population
            kids(a,dels(delsToRemove))=false;
        end
    end
end

end