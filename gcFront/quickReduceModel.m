function [model,delsToInds]=quickReduceModel(model, biomassRxn, targetRxn, minGrowth, minProd, tol, deleteGenes, delOptions)
% reduce the size of the input model by removing dead reactions, pooling unbranched
% pathways, and find the indices of a list of potential gene/reaction
% knockouts


m=convertModel(model);
m.lb(m.c==1)=minGrowth;
environment=getEnvironment();

if delOptions.skipReduction==0
    
    % finding reactions that cannot carry flux
    
    % first, solve model and identify any reactions that carry flux in the
    % solution
    s=optimizeCbModel2(m);
    if s.stat~=1
        error('Model is infeasible- try adjusting constraints or minimum growth value')
    else
        active=abs(s.x)>tol;
    end
    
    % test each reaction that wasn't carrying flux in the initial solution to
    % see if it can carry flux under any conditions
    
    parfor a=1:length(active)
        if active(a)
            continue
        end
        
        restoreEnvironment(environment,0);
        
        mTemp=m;
        mTemp.c=mTemp.c*0;
        mTemp.c(a)=1;
        
        s1=optimizeCbModel2(mTemp);
        if abs(s1.f)>tol
            active(a)=true;
            continue
        end
        
        mTemp.c(a)=-1;
        s2=optimizeCbModel2(mTemp);
        if abs(s2.f)>tol
            active(a)=true;
        end
    end
    
    % remove reactions that cant carry flux from the model
    model=removeRxns(model,model.rxns(~active));
    s2=optimizeCbModel(model);
    
    % check if removing the zero-flux reactions changed growth rate- if it
    % did, warn the user in case the reactions being removed actually do matter
    if abs(s2.f-s.f)>tol || isnan(s2.f)
        warning('on','all');
        warning(['Model reduction made growth rate shift by ', num2str(abs(s2.f-s.f)),' when minimum flux tolerance is ',num2str(tol)])
    end
    
    
    
    
    % check that target reaction wasn't removed from model
    if ~ismember(model.rxns,targetRxn)
        error('Target reaction cannot carry flux- try adjusting constraints, minimum growth or minimum flux tolerance')
    elseif ~ismember(model.rxns,biomassRxn)
        error('Biomass reaction cannot carry flux- try adjusting constraints or minimum flux tolerance')
    else
        biomassInd=find(ismember(model.rxns,biomassRxn));
        targetInd=find(ismember(model.rxns,targetRxn));
    end
    
    
    
    
    % now find reactions that form unbranched pathways and pool them
    
    a=1;
    changedAnything=0;
    while 1
        
        % check if reactions associated with a metabolite form an unbranched chain
        if sum(model.S(a,:)~=0)==2
            
            
            
            % find the reactions that can produce or consume the metabolite
            producingRxns=find(or( and(model.S(a,:)>0, model.ub'>0), and(model.S(a,:)<0, model.lb'<0) ));
            consumingRxns=find(or( and(model.S(a,:)>0, model.lb'<0), and(model.S(a,:)<0, model.ub'>0) ));
            
            % skip if one of the reactions is the biomass or target reaction
            if sum(ismember([producingRxns,consumingRxns],[biomassInd,targetInd]))~=0
                
                a=a+1;
                if a>length(model.mets)
                    if changedAnything==true
                        % repeat the process in case combining reactions made it
                        % possible to combine more reactions
                        
                        changedAnything=false;
                        a=1;
                    else
                        % no more reactions can be combined
                        break
                    end
                end
                continue
            end
            
            
            
            % if one reaction is reversible, then determine which
            % direction it is acting in. If both reactions are reversible, then
            % just pick one
            if length(producingRxns)==1 && length(consumingRxns)==2
                % if only one reaction can produce the metabolite, then this
                % reaction cannot consume the metabolite
                consumingRxns=consumingRxns(~ismember(consumingRxns,producingRxns));
            elseif length(consumingRxns)==1 && length(producingRxns)==2
                % if only one reaction can consume the metabolite, then this
                % reaction cannot produce the metabolite
                producingRxns=producingRxns(~ismember(producingRxns,consumingRxns));
            elseif length(producingRxns)==2 && length(consumingRxns)==2
                % the metabolite can be produced and consumed by either reaction,
                % so it is not possible to tell which direction they act in from
                % just their bounds and stoichiometry- just pick one at random to
                % be the producer
                if model.S(a,producingRxns(1))>0
                    producingRxns=producingRxns(1);
                    consumingRxns=consumingRxns(2);
                else
                    producingRxns=producingRxns(2);
                    consumingRxns=consumingRxns(1);
                end
            end
            
            
            % reverse direction of producing/consuming reaction if necessary
            if model.S(a,producingRxns)>0
                % reaction is acting in fwd direction
                prodStoich=model.S(:,producingRxns);
                prodBounds=[model.lb(producingRxns),model.ub(producingRxns)];
                prodC=model.c(producingRxns);
            else
                % reaction is in reverse order- rewrite reaction order
                prodStoich=model.S(:,producingRxns)*-1;
                prodBounds=[-model.ub(producingRxns),-model.lb(producingRxns)];
                prodC=-model.c(producingRxns);
            end
            if model.S(a,consumingRxns)<0
                % reaction is acting in fwd direction
                consStoich=model.S(:,consumingRxns);
                consBounds=[model.lb(consumingRxns),model.ub(consumingRxns)];
                consC=model.c(consumingRxns);
            else
                % reaction is in reverse order- rewrite reaction order
                consStoich=model.S(:,consumingRxns)*-1;
                consBounds=[-model.ub(consumingRxns),-model.lb(consumingRxns)];
                consC=-model.c(consumingRxns);
            end
            
            % adjust stoichiometry of consuming reaction to match producing
            % reaction
            multiplier=-prodStoich(a)/consStoich(a);
            consStoich=consStoich*multiplier;
            consBounds=consBounds/multiplier;
            consC=consC*multiplier;
            
            
            % create pooled reaction by combining stoichiometry + bounds
            jointStoich=prodStoich+consStoich;
            jointBounds=[max(prodBounds(1),consBounds(1)), min(prodBounds(2),consBounds(2))];
            jointC=prodC+consC;
            
            if jointBounds(2)==0
                % reaction is only going in reverse direction- reorder reaction
                jointStoich=-jointStoich;
                jointBounds=[-jointBounds(2),-jointBounds(1)];
                jointC=-jointC;
            end
            
            %          surfNet(model,model.mets(a))
            
            model.S(:,producingRxns)=jointStoich;
            model.S(:,consumingRxns)=0;
            model.lb(producingRxns)=jointBounds(1);
            model.ub(producingRxns)=jointBounds(2);
            model.c(producingRxns)=jointC;
            
            % combine fields associated with reactions
            
            % pool rxns
            
            %         disp(model.rxns([consumingRxns,producingRxns]))
            model.rxns{producingRxns}=[model.rxns{producingRxns},'/',model.rxns{consumingRxns}];
            %         disp(model.rxns(producingRxns))
            
            % pool rules
            
            
            
            if isfield(model,'rules')
                %                  disp(model.rules([consumingRxns,producingRxns]))
                if isempty(model.rules{consumingRxns})
                    % use rule from producing reaction
                elseif isempty(model.rules{producingRxns})
                    % use rule from consuming reaction
                    model.rules{producingRxns}=model.rules{consumingRxns};
                elseif ~isequal(model.rules{producingRxns},model.rules{consumingRxns})
                    % both reactions have distinct genes associated- combine the rules with an and statement
                    
                    if contains(model.rules{producingRxns},'|')
                        if contains(model.rules{consumingRxns},'|')
                            % put brackets around both rules- since combining
                            % rules with an and statement may otherwise lead to
                            % ambiguity
                            
                            model.rules{producingRxns}=['( ',model.rules{producingRxns},' ) & ( ', model.rules{consumingRxns}, ' )'];
                        else
                            % put brackets around only first rule
                            
                            model.rules{producingRxns}=['( ',model.rules{producingRxns},' ) & ', model.rules{consumingRxns}];
                        end
                    elseif contains(model.rules{consumingRxns},'|')
                        % put brackets only around second rule
                        
                        model.rules{producingRxns}=[model.rules{producingRxns},' & ( ', model.rules{consumingRxns}, ' )'];
                    else
                        % no brackets necessary- no ambiguity in rules
                        
                        model.rules{producingRxns}=[model.rules{producingRxns},' & ', model.rules{consumingRxns}];
                        
                    end
                    
                    % if both rules only contain and statements, then add rules
                    % together without a bracket
                    
                end
                
                %
                %              disp(model.rules(producingRxns))
                
            end
            
            % pool grRules
            if isfield(model,'grRules')
                %             disp(model.grRules([consumingRxns,producingRxns]))
                if isempty(model.grRules{consumingRxns})
                    % use rule from producing reaction
                elseif isempty(model.grRules{producingRxns})
                    % use rule from consuming reaction
                    model.grRules{producingRxns}=model.grRules{consumingRxns};
                elseif ~isequal(model.grRules{producingRxns},model.grRules{consumingRxns})
                    % both reactions have distinct genes associated- combine the rules with an and statement
                    model.grRules{producingRxns}=['( ',model.grRules{producingRxns},' ) and ( ', model.grRules{consumingRxns}, ' )'];
                end
                
                
                %             disp(model.grRules(producingRxns))
            end
            
            % pool rxnGeneMat
            if isfield(model,'rxnGeneMat')
                %             disp(find(model.rxnGeneMat(consumingRxns,:)))
                %             disp(find(model.rxnGeneMat(producingRxns,:)))
                
                genesInCons=model.rxnGeneMat(consumingRxns,:)~=0;
                model.rxnGeneMat(producingRxns,genesInCons)=1;
                
                %             disp(find(model.rxnGeneMat(producingRxns,:)))
            end
            
            % pool rxnNames
            if isfield(model,'rxnNames')
                
                %             disp(model.rxnNames([consumingRxns,producingRxns]))
                if isempty(model.rxnNames{consumingRxns})
                    % use name from producing reaction
                elseif isempty(model.rxnNames{producingRxns})
                    % use name from consuming reaction
                    model.rxnNames{producingRxns}=model.rxnNames{consumingRxns};
                else
                    % both reactions have names associated- combine the names
                    % with a slash
                    model.rxnNames{producingRxns}=[model.rxnNames{producingRxns}, '/', model.rxnNames{consumingRxns}];
                end
                
                %             disp(model.rxnNames(producingRxns))
            end
            
            % pool subSystems
            if isfield(model,'subSystems')
                %             disp(model.subSystems{consumingRxns})
                %             disp(model.subSystems{producingRxns})
                if isempty(model.subSystems{consumingRxns})
                    % use subsystem from producing reaction
                elseif isempty(model.subSystems{producingRxns})
                    % use subsystem from consuming reaction
                    model.subSystems{producingRxns}=model.subSystems{consumingRxns};
                elseif ~isequal(model.subSystems{producingRxns},model.subSystems{consumingRxns})
                    % both reactions have distinct subsystems associated- combine the
                    % subsystems with a slash
                    model.subSystems{producingRxns}=join([model.subSystems{producingRxns}, '/', model.subSystems{consumingRxns}],'');
                end
                
                
                %             disp(model.subSystems{producingRxns})
            end
            
            
            
            %          surfNet(model,model.rxns(producingRxns))
            %         disp(' ')
            
            % something was changed during this iteration of the loop
            changedAnything=true;
        end
        a=a+1;
        if a>length(model.mets)
            if changedAnything==true
                % repeat the process in case combining reactions made it
                % possible to combine more reactions
                
                changedAnything=false;
                a=1;
            else
                % no more reactions can be combined
                break
            end
        end
        
    end
    
    
    deletedRxns=sum(abs(model.S),1)==0;
    model=removeRxns(model,model.rxns(deletedRxns));
    
    
    s3=optimizeCbModel(model);
    % check if pooling reactions affected growth rate at all
    if abs(s3.f-s.f)>tol || isnan(s3.f)
        warning(['Reaction pooling during model reduction made growth rate shift by ', num2str(abs(s3.f-s.f)),' when tolerance to flux deviation has been set to ',num2str(tol),' - something may have gone wrong during model reduction.'])
    end
    
end



% if deleting reactions, then test reactions for if they are appropriate
% knockouts
if ~deleteGenes
    % create a lpModel for use in parfor loop
    lpModel=convertModel(model);
%     set minimum bound on product so FBA will only have a non-zero
%     objective if KO does not prevent product synthesis
   lpModel.lb(ismember(model.rxns,targetRxn))=minProd;

%     % set minimum bound on growth and change objective to product
%     % synthesis, so FBA will only have a non-zero objective if KO does not 
%     % prevent growth or product synthesis
%     lpModel.lb(lpModel.c==1)=minGrowth;
%     lpModel.c(ismember(model.rxns,biomassRxn))=0;
%     lpModel.c(ismember(model.rxns,targetRxn))=1;
    
    % create vector of valid reaction deletions
    validDels=~ismember(model.rxns,{targetRxn,biomassRxn});
    
    % exclude reactions that are not gene associated
    if delOptions.onlyKoGeneAssoc
        validDels(ismember(model.rules,''))=false;
    end
    
    % exclude reactions that must have flux
    validDels(model.lb>0)=false;
    validDels(model.ub<0)=false;
    
    
    % exclude reactions that are on the list of reactions to be ignored
    if ~isempty(delOptions.ignoreListRxns)
        for a=1:length(model.rxns)
            if validDels(a)==true
                
                % check if reaction is part of the list of reactions to be
                % ignored
                isForbidden=ismember(splitString(model.rxns{a},'/'), delOptions.ignoreListRxns);
                % exclude if it is on list
                if sum(isForbidden)==length(isForbidden)
                    validDels(a)=false;
                end
            end
        end
    end
    
    
    % exclude reactions if they cannot be knocked out without knocking out
    % an essential gene or a gene on the list of genes to not be knocked
    % out
    if delOptions.dontKoEss || ~isempty(delOptions.ignoreListGenes)
        % define excluded reactions as essential, so that their deletion is
        % not permissible
        if ~isempty(delOptions.ignoreListGenes)
            essGeneList=ismember(model.genes,delOptions.ignoreListGenes);
        else
            essGeneList=false(size(model.genes));
        end
        
        if delOptions.dontKoEss
            % find computationally essential genes and exclude any reaction
            % that can't be knocked out without KO of a computationally
            % essential gene
            koInds=false(length(model.genes),length(model.rxns));
            % establish which knockouts occur after gene KO
            for a=1:length(essGeneList)
                if essGeneList(a)==true
                    % skip if we are already not considering this gene
                    continue
                end
                
                x=true(size(model.genes));
                x(a)=false;
                
                delInds=find(model.rxnGeneMat(:,a)~=0);
                
                for b=1:length(delInds)
                    koInds(a,delInds(b))=~eval(model.rules{delInds(b)});
                end
                %                 disp(model.genes(a))
                %                 disp(model.grRules(logical(model.rxnGeneMat(:,a))))
                %                 disp(koInds(a,delInds))
                %                 disp(' ')
            end
            
            
            % remove non-unique KO combinations so same optimisation isnt
            % carried out twice
            
            [koInds,~,iC]=unique(koInds,'rows');
            invalidDeletion=false(size(koInds,1),1);
            
            parfor a=1:size(koInds,1)
                
                restoreEnvironment(environment,0);
                
                
                m=lpModel;
                m.lb(koInds(a,:))=0;
                m.ub(koInds(a,:))=0;
                sTemp=optimizeCbModel2(m);
                if ~(minGrowth<=sTemp.f)
                %if ~(minProd<=sTemp.f)
                    % model is infeasible or growth/product are below
                    % threshold
                    invalidDeletion(a)=true;
                    
                end
                
            end
            
            essGeneList=or(essGeneList,invalidDeletion(iC));
            
        end
        
        
        % if a single reaction KO was found to not be feasible, then remove the
        % reaction from the list of reactions to be considered (saves
        % calculating this again later)
        singleInfeasible=and(invalidDeletion,sum(koInds,2)==1);
        [~,lethalDel]=find(koInds(singleInfeasible,:));
        validDels(lethalDel)=false;
        
        % test if reaction remains active if all non-essential genes are
        % knocked out. If it does, then it can't be knocked out without an
        % essential gene ko, and so shouldnt be considered further
        x=essGeneList;
        for a=1:length(model.rxns)
            if validDels(a) && ~isempty(model.rules{a})
                validDels(a)=~eval(model.rules{a});
            end
        end
    end
    
    
    % exclude deletions if they prevent growth above threshold/product
    % synthesis
    parfor a=1:length(model.rxns)
        if ~validDels(a)
            % reaction is not gene associated, or reaction is
            % growth/product synthesis
            continue
        end
        
        restoreEnvironment(environment,0);
        
        % test deletion of reaction
        m=lpModel;
        m.lb(a)=0;
        m.ub(a)=0;
        sTemp=optimizeCbModel2(m);
        
        if ~(minGrowth<=sTemp.f)
        %if ~(minProd<=sTemp.f)
            % model is infeasible or growth/product is below threshold
            validDels(a)=false;
        end
        
    end
    
    delsToInds=find(validDels);
    
    
    
else
    % deleting genes
    
    % Find genes with no associated reaction
    unassocGenes=find(sum(full(model.rxnGeneMat),1)==0);
    
    % update the numbering of the rules field for when these genes are
    % removed
    counter=0;
    for a=1:length(model.genes)
        
        if ismember(a,unassocGenes)
            % gene has been removed- update counter
            counter=counter+1;
        elseif counter~=0
            % the number of the gene has changed, so rules field must be updated
            assocRxns=find(full(model.rxnGeneMat(:,a))~=0);
            for b=1:length(assocRxns)
                model.rules{assocRxns(b)}=strrep(model.rules{assocRxns(b)},['x(',num2str(a),')'],['x(',num2str(a-counter),')']);
            end
        end
    end
 
    % delete the genes and their field in rxngenemat
    model.genes(unassocGenes)=[];
    model.rxnGeneMat(:,unassocGenes)=[];

    
    % Pooling genes
    
    geneMat=full(model.rxnGeneMat);
    
    [~,~,iC]=unique(geneMat','rows');
    
    combined=zeros(length(model.genes),1);
    
    % Combine genes that form equivalent sets where any one could be
    % deleted for the same effect
    
    % loop over genes
    for a=1:length(model.genes)
        
        
        % skip if this gene has been combined already
        if combined(a)~=0
            continue
        end
        
        % find other genes that have matching associated rxns and havent already been
        % combined
        matchingInds=setdiff( find( and( combined==0, iC==iC(a) ) ), a);
        
        % stop if there are no genes with matching reactions
        if isempty(matchingInds)
            continue
        end

        % find the index of the reactions associated with these genes
        assocRxnInds=find(geneMat(:,a));
       
        orInds=true(size(matchingInds));
        
        % loop over rxns associated with our gene
        for b=1:length(assocRxnInds)
            
            ruleMat=readRulesField(model.rules{assocRxnInds(b)});
           
            % find columns which the initial gene is in
            currInds=logical(sum(ruleMat==a,1));
          
            for c=1:length(matchingInds)
                % find columns which the other genes are in
                tempInds=logical(sum(ruleMat==matchingInds(c),1));
                
                % if columns don't match, then deletion of these genes
                % can have different effects, so their deletions are
                % not equivalent. Deletions will only be equivalent if we
                % have identical columns for every reaction
                if ~isequal(currInds,tempInds)
                    orInds(c)=false;
                end
                
            end
            
        end
        
        % combining genes into equivalent deletion sets
        model.genes{a}=char(join(sort([model.genes(a);model.genes(matchingInds(orInds))]),'/'));
        
        % mark combined genes so they aren't analysed further
        combined(matchingInds(orInds))=1;
        
    end
    
    % find sets where all of the genes must be deleted to have any effect
    
    % loop over genes
    for a=1:length(model.genes)
        
        % skip if this gene has been combined already
        if combined(a)~=0
            continue
        end
        
        % find other genes that have matching associated rxns and havent already been
        % combined
        matchingInds=setdiff( find( and( combined==0, iC==iC(a) ) ), a);
        
        % stop if there are no genes with matching reactions
        if isempty(matchingInds)
            continue
        end
        
        % find the index of the reactions associated with these genes
        assocRxnInds=find(geneMat(:,a));
        
        andInds=true(size(matchingInds));
        
        % loop over rxns associated with our gene
        for b=1:length(assocRxnInds)
            ruleMat=readRulesField(model.rules{assocRxnInds(b)});
            
            % find columns which the initial gene is in
            currInds=logical(sum(ruleMat==a,1));
            
            for c=1:length(matchingInds)
                % find columns which the other genes are in
                tempInds=logical(sum(ruleMat==matchingInds(c),1));
                
                % if columns match, then deletion of these genes has the
                % same effect, so we should not be deleting them
                % simultaneously
                if isequal(currInds,tempInds)
                    andInds(c)=false;
                else
                    
                    % we should also check to see if these genes are part
                    % of identical sets, since if they are they should be
                    % knocked out simultaneously
                    
                    % get gene matrix, and replace initial gene number with
                    % current gene number
                    aaa=ruleMat;
                    aaa(ruleMat==a)=matchingInds(c);
                    
                    % get columns that correspond to initial and current
                    % gene
                    bbb=aaa(:,currInds);
                    ccc=aaa(:,tempInds);
                    
                    % sort columns by number
                    bbb=sortrows(bbb')';
                    ccc=sortrows(ccc')';
                    
                    % check if the columns that initial and current gene
                    % are in are identical except for them- if they are
                    % not, we shouldn't knock out these genes
                    % simultaneously as they may have different effects
                    if ~isequal(bbb,ccc)
                        andInds(c)=false;
                    end
                end
                
            end
            
        end
        
        % combining genes into equivalent deletion sets
        model.genes{a}=char(join(sort([model.genes(a);model.genes(matchingInds(andInds))]),'+'));
        
        % mark combined genes so they aren't analysed further
        combined(matchingInds(andInds))=2;
        
        
    end
    
    % remove all genes that have been combined, and rewrite rules and
    % grRules fields
    
    orCombinedInds=find(combined==1);
    andCombinedInds=find(combined==2);
    
    % track the genes that are not actually contributing to a reaction
    unusedGenes=true(size(model.genes));
    
    for a=1:length(model.rxns)
        
        
        testRule=model.rules{a};
        
        if isempty(testRule)
            finalRule='';
        else
            % get the old rule and convert it into a matrix of gene
            % indices, where columns represent combinations of genes that
            % will sustain a reaction
            ruleMat=readRulesField(testRule);
            
            % remove genes that have already been combined from ruleMat- so they won't
            % be included in the newly written rule
            [~,yInd]=find(ismember(ruleMat,andCombinedInds));
            ruleMat(:,yInd)=[];
            ruleMat(ismember(ruleMat,orCombinedInds))=inf;
            
            % remove columns of grRules if they contain all the genes in
            % another column, since this column is redundant. For example,
            % if grRule is: A or A+B or A+C, then A+B and A+C are
            % redundant when considering KOs, as the only way the reaction
            % can be KOd is if A is knocked out. Thus, there is no point in
            % retaining this part of the grRule for our reduced model
            if size(ruleMat,2)>1
                % create list of columns to be deleted
                columnsToDel=false(1,size(ruleMat,2));
                
                for b=1:size(ruleMat,2)
                    
                    % don't analyse column if it is set to be deleted
                    if columnsToDel(b)==true
                        continue
                    end
                    
                    % create temporary list of columns that contain same
                    % genes as current column
                    foundCols=true(size(columnsToDel));
                    foundCols(b)=false;
                    
                    % find genes in current column
                    columnIndices=setdiff(unique(ruleMat(:,b)),inf);
                    
                    % check other columns to find ones that contain
                    % same genes as current column
                    for c=1:length(columnIndices)
                        foundCols=and(foundCols,sum(ruleMat==columnIndices(c),1)~=0);
                        if sum(foundCols)==0
                            break
                        end
                    end
                    
                    % add any column that contained the same genes to the
                    % list of columns to delete
                    columnsToDel=or(columnsToDel,foundCols);
                end
                
                % delete all redundant columns
                ruleMat=ruleMat(:,~columnsToDel);

            end
            
            
            % construct new grRules from the rule matrix
            
            
            % find genes found in all sets
            uniqueGenes=setdiff(unique(ruleMat),inf);
            
            % record that these genes are being used so they arent deleted later
            unusedGenes(uniqueGenes)=false;
            
            % from list of genes, find genes that are present in every
            % column of rule matrix
            inAllSets=false(size(uniqueGenes));
            for b=1:length(uniqueGenes)
                % find columns that the gene is in
                uniqueGenesCount=sum(ruleMat==uniqueGenes(b),1);
                % if it is in all columns, then add it to inAllSets
                if min(uniqueGenesCount)>=1
                    inAllSets(b)=true;
                end
            end
            
            % make a string of the genes in every set, joined by AND
            % statements
            if sum(inAllSets)~=0
                ruleString1=char(join(model.genes(uniqueGenes(inAllSets)),' and '));
            else
                ruleString1='';
            end
            
            % now make a string of everything else
            if sum(inAllSets)~=length(inAllSets)
                
                % remove anything that has already been written into the gene rule
                ruleMat(ismember(ruleMat,uniqueGenes(inAllSets)))=inf;
                
                % remove any columns that are now non-unique
                ruleMat=sort(ruleMat);
                ruleMat=unique(ruleMat','rows')';
                
                % remove any column with no gene index in it
                onlyInf=sum(ruleMat==inf,1)==size(ruleMat,1);
                ruleMat=ruleMat(:,~onlyInf);
                
                
                % join together the genes in each column of ruleMat with AND statements
                tempRules=strings(1,size(ruleMat,2));
                for b=1:size(ruleMat,2)
                    tempRow=unique(ruleMat(:,b));
                    tempRow=setdiff(tempRow,inf);
                    if ~isempty(tempRow)
                        tempRules(b)=join(model.genes(tempRow),' and ');
                    end
                end
                
                
                % now join together the different columns with OR statements
                
                % join together the columns with only 1 gene (don't add a bracket)
                onlyOneGene=sum(ruleMat~=inf,1)==1;
                if sum(onlyOneGene)~=0
                    tempString1=char(join(tempRules(onlyOneGene),' or '));
                else
                    tempString1='';
                end
                % join together the columns with more than one gene (add brackets)
                if sum(onlyOneGene)~=length(onlyOneGene)
                    tempString2=['( ',char(join(tempRules(~onlyOneGene),' ) or ( ')),' )'];
                else
                    tempString2='';
                end
                
                % now combine these strings (with OR statement if
                % necessary)
                if isempty(tempString1)
                    ruleString2=tempString2;
                elseif isempty(tempString2)
                    ruleString2=tempString1;
                else
                    ruleString2=[tempString1,' or ', tempString2];
                end
                
                
            else
                ruleString2='';
            end
            
            % now combine the genes that were in everything with the genes
            % that only appeared in some columns
            
            if isempty(ruleString1)
                % no need to join the first string to second if there is no first
                % string
                finalRule=ruleString2;
            elseif isempty(ruleString2)
                % no need to join the second string to first if there
                % is no second string
                finalRule=ruleString1;
            else
                finalRule=[ruleString1,' and ( ',ruleString2,' )'];
            end
        end
        
        % store the rewritten rule in the grRules field
        model.grRules{a}=finalRule;
        
    end

    
    % remove genes that were combined into another gene, or that do not
    % contribute to a reaction after redundancy was removed from this
    % reaction
    model.genes=model.genes(and(combined==0,~unusedGenes));
    model.rxnGeneMat=model.rxnGeneMat(:,and(combined==0,~unusedGenes));
    
    % rewrite rules field based on the grRules field
    model.rules=strrep(model.grRules,'or','|');
    model.rules=strrep(model.rules,'and','&');
    for a=1:length(model.genes)
        model.rules=strrep(model.rules,model.genes{a},['x(',num2str(a),')']);
    end

    
    % create list of genes
    validDels=true(size(model.genes));
    
    
    % remove genes linked to reactions to be ignored, and reactions that
    % must have flux
    
    for a=1:length(model.rxns)
        
        % don't allow deletion of genes linked to reactions that must have
        % flux, as deletion of these genes may cause problems later. Could
        % probably make this work by having a list of all reactions that
        % must carry flux, and then not allowing gene knockout combinations
        % that result in KO of a reaction that must have flux.
        if model.lb(a)>0 || model.ub(a)<0
            linkedGenes=logical(full(model.rxnGeneMat(a,:)));
            validDels(linkedGenes)=false;
            continue
        end
        
        if ~isempty(delOptions.ignoreListRxns)
            isForbiddenRxn=ismember(splitString(model.rxns{a},'/'), delOptions.ignoreListRxns);
            if sum(isForbiddenRxn)==length(isForbiddenRxn)
                linkedGenes=logical(full(model.rxnGeneMat(a,:)));
                validDels(linkedGenes)=false;
            end
        end
        
    end
    
    
    % check for genes on the list of genes to be ignored
    if ~isempty(delOptions.ignoreListGenes)
        for a=1:length(validDels)
            % check if gene is part of the list of reactions to be
            % ignored
            isForbidden=ismember(splitString(model.genes{a},'/'), delOptions.ignoreListGenes);
            % exclude if it is on list
            if sum(isForbidden)==length(isForbidden)
                validDels(a)=false;
            end
        end
    end
    
    
    % create matrix of reaction deletions
    rxnDelMat=false(length(validDels),length(model.rxns));
    for a=1:length(validDels)
        if validDels(a)==true
            assocRxns=find(full(model.rxnGeneMat(:,a)));
            if length(assocRxns)==1
                % deletion of gene will either do nothing, or KO the only
                % assoc rxn, so may as well test if this rxn is essential
                rxnDelMat(a,assocRxns)=true;
            else
                % create variable x for evaluating the rules statements
                x=true(length(model.genes),1);
                % KO the gene of interest in x
                x(a)=false;
                % test each associated reaction to see if it is KOd if the
                % gene of interest is KOd, and store result in rxnDelMat
                for b=1:length(assocRxns)
                    rxnDelMat(a,assocRxns(b))=~eval(model.rules{assocRxns(b)});
                end
            end
        end
    end
    
    [rxnDelMat,~,newToOld]=unique(rxnDelMat,'rows');
    feasibleSolutions=true(size(rxnDelMat,1),1);
    
    % create a lpModel for use in parfor loop
    lpModel=convertModel(model);
    lpModel.lb(ismember(model.rxns,targetRxn))=minProd;
%     % set minimum bound on growth and change objective to  product synthesis
%     lpModel.lb(lpModel.c==1)=minGrowth;
%     lpModel.c(ismember(model.rxns,biomassRxn))=0;
%     lpModel.c(ismember(model.rxns,targetRxn))=1;
    
    % find if the reaction KOs caused by gene KOs stop model from meeting
    % growth/product thresholds
    parfor a=1:size(rxnDelMat,1)
        restoreEnvironment(environment,0);
        % test gene knockouts to see if they prevent growth/product
        % synthesis
        
        m=lpModel;
        m.lb(rxnDelMat(a,:))=0;
        m.ub(rxnDelMat(a,:))=0;
        sTemp=optimizeCbModel2(m);
        
        if ~(minGrowth<=sTemp.f)
            feasibleSolutions(a)=false;
        end
    end
    
    % remove any gene deletion if the reaction deletions it is linked to
    % lower growth/product by too much
    feasSols2=feasibleSolutions(newToOld);

    validDels=and(validDels,feasSols2);
    
    
    delsToInds=find(validDels);
    
    
end


end