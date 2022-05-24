function kids=myCrossover(parents, ~, genomeLength, ~, ~, population, maxKnockouts, geneCounter)    
% carry out crossovers during multiobjective GA

kids=false(length(parents)/2,genomeLength);
for a=1:size(kids,1)
    for b=1:genomeLength
        if rand>0.5
            kids(a,b)=population( parents(a*2-1), b );
        else
            kids(a,b)=population( parents(a*2), b );
        end
    end
end

% check if child has too many knockouts

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