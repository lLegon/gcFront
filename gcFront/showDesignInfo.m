function txt=showDesignInfo(~, info, designTable)
% display the metrics of a design that the user clicked on in interactive 
% Pareto, and make the "prod env" and "show secretion" buttons link to this design

% check what user clicked on
if isa(info.Target,'matlab.graphics.chart.primitive.Scatter')
    
    
    % they clicked on a datapoint- show metrics of point
    if max(designTable{:,4})>0
        % coupled designs have been found
        x = info.Position(1);
        y = info.Position(2);
        
        ind=find(and(designTable{:,3}==x, designTable{:,4}==y));
    else
        x=info.Position(1);
        y='0';
        
        ind=find(and(designTable{:,3}==x, designTable{:,5}==info.Position(2)));
    end
    
    % find coupling strength for design(s)
    couplingNum=designTable{ind,5};
    if length(unique(couplingNum))~=1
        % there are multiple designs with this combination of
        % growth/product synthesis and they have unique coupling strength
        % values
        couplingVal=join(string(couplingNum),' / ');
    else
        couplingVal=string(couplingNum(1));
    end
    
    % find number of deletions for design(s)
    delNum=designTable{ind,2};
    if length(unique(delNum))~=1
        % there are multiple designs with this combination of
        % growth/product synthesis and they have unique numbers of
        % deletions
        delVal=join(string(delNum),' / ');
    else
        delVal=string(delNum(1));
    end
    
    % finding distance of design to ideal point in obj space
    dist = designTable{ind,6};
    if length(unique(dist)) ~= 1
        distVal = join(string(dist),' / ');
    else
        distVal = string(dist(1));
    end
    
    if ~isequal(y,'0')
        txt=["Deletions:";
            designTable{ind,1};
            join(["Number of deletions = ", delVal],'');
            join(["Growth = ", x],'');
            join(["Product synthesis = ", y],'');
            join(["Coupling strength = ", couplingVal],'');
            join(["Distance to ideal = ", distVal],'')];
    else
        txt=["Deletions:";
            designTable{ind,1};
            join(["Number of deletions = ", delVal],'');
            join(["Growth = ", x],'');
            join(["Product synthesis = ", y],'');
            join(["Coupling strength = ", couplingVal],'')];
    end
    
    % find and put deletions into the callback for the production envelope and show secretion buttons
    currObj=findobj(gcf,'type','UIControl');
    
    for a=1:length(currObj)
        if strcmp(currObj(a).String,'Prod Env') || strcmp(currObj(a).String,'Show secretion')
            currObj(a).Callback{3}=designTable{ind,1};
        end
    end
    
    
    
    
elseif isa(info.Target,'matlab.graphics.chart.primitive.Line')
    % they clicked on the production envelope- show max growth and product
    % synthesis
    txt=["WT Production Envelope";
        join(["Max growth = ",num2str(max(info.Target.XData))],'');
        join(["Max product synthesis = ",num2str(max(info.Target.YData))],'')];
    
    % put a blank into showSecretion and ProdEnv
    currObj=findobj(gcf,'type','UIControl');
    
    for a=1:length(currObj)
        if strcmp(currObj(a).String,'Prod Env') || strcmp(currObj(a).String,'Show secretion')
            currObj(a).Callback{3}='';
        end
    end
end


end