function prodEnvButtonFunction(~,~,model,deletions,targetRxn, wtEnv, deleteGenes)
% When user presses "prodEnv" button, display a production envelope for the
% most recently selected design

if ~isempty(deletions)
    while 1
        % get current figures
        currFigures=findobj('type','figure');
        
        % get index of the current figures
        %paretoInd=[];
        prodEnvInd=[];
        
        for a=1:length(currFigures)
            if strcmp(currFigures(a).Name,'Production Envelopes')
                prodEnvInd=a;
%             elseif strcmp(currFigures(a).Name,'Pareto Front')
%                 paretoInd=a;
            end
        end
        
        if isempty(prodEnvInd)
            % Create a figure for production envelopes if one does not already
            % exist
            figure('Name','Production Envelopes')
            plot(wtEnv(:,1),wtEnv(:,2),'k','LineWidth',3)
            legend('WT')
            xlabel('Growth rate (/h)')
            ylabel('Target flux (mmol/gDW/h)')
            ylim([0,ceil(max(wtEnv(:,2)))])
            box off
        else
            % a production envelope figure exists- continue
            break
        end
    end
    
    figure(currFigures(prodEnvInd))
    
    ax=gca;
    if isprop(ax.Legend,'String')
        legendActive=true;
        origLegend=ax.Legend.String;
        
        % skip deletions that have already been plotted
        deletions=deletions(~ismember(deletions,origLegend));
    else
        legendActive=false;
    end
    
    % find and plot production envelope for each input deletion
    hold on
    for a=1:size(deletions,1)
    
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
        tempEnv=prodEnvFast(m,targetRxn,0);
        plot(tempEnv(:,1),tempEnv(:,2),'LineWidth',3, 'DisplayName', deletions(a,:))
    end
    hold off
    
    if legendActive
        legend([origLegend,deletions'],'Interpreter','none');
    end

    
    %figure(currFigures(paretoInd))
end

end