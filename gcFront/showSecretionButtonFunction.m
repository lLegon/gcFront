function showSecretionButtonFunction(~,~,model,deletions, deleteGenes, excRxns)
% when user presses showSecretion, show flux through exchange reactions 
% after deletion of specified reactions


if ~isempty(deletions)
    for a=1:size(deletions,1)
        
        disp(['Deletions: ',char(deletions(a))])
        
        m=model;
        if ~ismember(deletions{a},'')
            if ~deleteGenes
                % delete reactions corresponding to deletions
                dels=splitString(deletions{a},' ');
                dels2=splitString(dels,'/');
                for b=1:length(dels2)
                    dels2{b}=dels2{b}{1};
                end
                
                m.lb(ismember(m.rxns,dels2))=0;
                m.ub(ismember(m.rxns,dels2))=0;
                
            else
                % find genes corresponding to deletions
                delGenes=splitString(deletions{a},' ');
                delGenes2=splitString(delGenes,'/');
                for b=1:length(delGenes2)
                    delGenes2{b}=delGenes2{b}{1};
                end
                
                x=~ismember(m.genes,delGenes2);
                
                % delete reactions knocked out by these gene deletions
                assocRxns=find(sum(m.rxnGeneMat(:,~x),2)>0);
                delRxnInds=false(size(model.rxns));
                for b=1:length(assocRxns)
                    delRxnInds(assocRxns(b))=~eval(model.rules{assocRxns(b)});
                end
                m.lb(delRxnInds)=0;
                m.ub(delRxnInds)=0;
                
            end
        end
        
        
        if sum(excRxns)>0
            % find min/max flux of exchange reactions at max growth
            [lowers,uppers]=fluxVariability(m,100,'max',m.rxns(excRxns));
            
            % find metabolites linked to exchange reactions
            excMetNames=strings(size(excRxns,1),1);
            for b=1:length(excRxns)
                excMet=model.S(:,excRxns(b))~=0;
                if isfield(model,'metNames')
                    excMetNames(b)=join(model.metNames(excMet),', ');
                else
                    excMetNames(b)=join(model.mets(excMet),', ');
                end
            end
            
            % find formulas of exchange reactions
            formulas=printRxnFormula(model,'rxnAbbrList',model.rxns(excRxns),'printFlag',false);
            
            
            
            t=table(excMetNames,lowers,uppers,model.rxns(excRxns),formulas,'VariableNames',{'Metabolite','MinFlux','MaxFlux','Reaction','ReactionFormula'});
            t=sortrows(t,3,'descend');
            active=or(t{:,2}~=0,t{:,3}~=0);
            disp(t(active,:))
        end
    end
else
    disp('Deletions:')
    
    % just find secretion of WT model
    
    % find min/max flux of exchange reactions at max growth
    [lowers,uppers]=fluxVariability(model,100,'max',model.rxns(excRxns));
    
    % find metabolites linked to exchange reactions
    excMetNames=strings(size(excRxns,1),1);
    for b=1:length(excRxns)
        excMet=model.S(:,excRxns(b))~=0;
        if isfield(model,'metNames')
            excMetNames(b)=join(model.metNames(excMet),', ');
        else
            excMetNames(b)=join(model.mets(excMet),', ');
        end
    end
    
    % find formulas of exchange reactions
    formulas=printRxnFormula(model,'rxnAbbrList',model.rxns(excRxns),'printFlag',false);
    
    
    
    t=table(excMetNames,lowers,uppers,model.rxns(excRxns),formulas,'VariableNames',{'Metabolite','MinFlux','MaxFlux','Reaction','ReactionFormula'});
    t=sortrows(t,3,'descend');
    active=or(t{:,2}~=0,t{:,3}~=0);
    disp(t(active,:))
end

end