function [ prodEnvVals ] = prodEnvFast( model, targetInd, drawOption, tol, nVals, lineTol)
%[ prodEnvVals ] = prodEnv( model, targetInd, drawOption )
% returns production envelope values for a given model
% can supply a lp model and the logical index of the target reaction, or
% can supply the full model and the name of the reaction

if isfield(model,'rxns')
    % have not been supplied with an lp model
    if ~islogical(targetInd) && ~isnumeric(targetInd)
        targetInd=ismember(model.rxns,targetInd);
    end
    model=convertModel(model);
end
if ~islogical(targetInd)
    targetInd=logical(targetInd);
end
if nargin<3 || isempty(drawOption)
    drawOption=1;
end
if nargin<4 || isempty(tol)
    tol=10^-6;
end
if nargin<5 || isempty(nVals)
    nVals=500;
end
if nargin<6 || isempty(lineTol)
    lineTol=10^-3;
end

% minimise growth
m=model;
m.c=-m.c;
s=optimizeCbModel2(m);
% objective values become negative when minimising with optimizeCbModel2-
% so take values from flux vector
if s.stat==1
    s.f=s.x(model.c==1);
end

% minimise product at minimum growth
m.c=-double(targetInd);
m.ub(model.c==1)=s.f+tol;
s2=optimizeCbModel2(m);
% objective values become negative when minimising with optimizeCbModel2- so take values from
% flux vector
if s2.stat==1
    s2.f=s2.x(targetInd);
else
    warning('Production envelope is inaccurate- model was infeasible when calculating minimum product at min growth. Consider increasing tolerance to numerical error.')
    s2.f=0;
end

% maximise growth
m=model;
s3=optimizeCbModel2(m);

% minimise product at maximum growth
m=model;
m.c=-double(targetInd);
m.lb(model.c==1)=s3.f-tol;
s4=optimizeCbModel2(m);
% objective values become negative when minimising with optimizeCbModel2- so take values from
% flux vector
if s4.stat==1
    s4.f=s4.x(targetInd);
else
    warning('Production envelope is inaccurate- model was infeasible when calculating minimum product at max growth. Consider increasing tolerance to numerical error.')
    s4.f=0;
end

valueListL=[s.f,s2.f,0;s3.f,s4.f,1];

% set up temporary model that will minimise target synthesis
m=model;
m.c=-double(targetInd);
while length(valueListL)<inf
    newValueList=nan(size(valueListL));
    for a=1:size(valueListL,1)-1
        if valueListL(a,3)==0
            % test if the production envelope is a straight line between
            % two points. If it is, then no need to test any more points
            % along this segment of the production envelope.
            midGrowth=mean(valueListL(a:a+1,1));
            m.lb(model.c==1)=midGrowth-tol;
            m.ub(model.c==1)=midGrowth+tol;
            s5=optimizeCbModel2(m);
            % objective values become negative when minimising with optimizeCbModel2- so take values from
            % flux vector
            if s5.stat==1
                s5.f=s5.x(targetInd);
            else
                % min product could not be calculated- will assume that
                % production envelope is a straight line
                warning('Production envelope may be inaccurate- model was infeasible when calculating minimum product. Consider increasing tolerance to numerical error.')
                s5.f=mean(valueListL(a:a+1,2));
            end
            
            if abs(s5.f-mean(valueListL(a:a+1,2)))<lineTol
                valueListL(a,3)=1;
            else
                newValueList(a,:)=[midGrowth,s5.f,0];
            end
        end
    end
    newValueList=rmmissing(newValueList);
    if isempty(newValueList)
        break
    else
        valueListL=[valueListL;newValueList];
        valueListL=sortrows(valueListL,1);
    end
    
    if size(valueListL,1)>nVals
        break
    end
end

% minimise growth - can use previous value
% s=optimizeCbModel(model,'min');

% maximise product at minimum growth
m=model;
m.c=double(targetInd);
m.ub(model.c==1)=s.f+tol;
s2=optimizeCbModel2(m);
if s2.stat~=1
    warning('Production envelope is inaccurate- model was infeasible when calculating maximum product at min growth. Consider increasing tolerance to numerical error.')
    s2.f=0;
end

% maximise growth- can use previous value
% s3=optimizeCbModel(model,'max');

% maximise product at maximum growth
m=model;
m.c=double(targetInd);
m.lb(model.c==1)=s3.f-tol;
s4=optimizeCbModel2(m);
if s4.stat~=1
    warning('Production envelope is inaccurate- model was infeasible when calculating maximum product at max growth. Consider increasing tolerance to numerical error.')
    s4.f=0;
end

valueListH=[s.f,s2.f,0;s3.f,s4.f,1];

% create temporary model that maximises product
m=model;
m.c=double(targetInd);
while length(valueListH)<inf
    newValueList=nan(size(valueListH));
    for a=1:size(valueListH,1)-1
        if valueListH(a,3)==0
            midGrowth=mean(valueListH(a:a+1,1));
            m.lb(model.c==1)=midGrowth-tol;
            m.ub(model.c==1)=midGrowth+tol;
            s5=optimizeCbModel2(m);
            
            if s5.stat~=1
                warning('Production envelope may be inaccurate- model was infeasible when calculating maximum product. Consider increasing tolerance to numerical error.')
                s5.f=mean(valueListH(a:a+1,2));
            end
            
            if abs(s5.f-mean(valueListH(a:a+1,2)))<lineTol
                valueListH(a,3)=1;
            else
                newValueList(a,:)=[midGrowth,s5.f,0];
            end
        end
    end
    newValueList=rmmissing(newValueList);
    if isempty(newValueList)
        break
    else
        valueListH=[valueListH;newValueList];
        valueListH=sortrows(valueListH,1);
    end
    
    if size(valueListH,1)>nVals
        break
    end
end

prodEnvVals=[valueListL(:,[1,2]);valueListH((end:-1:1),[1,2])];
prodEnvVals=[prodEnvVals;prodEnvVals(1,:)]; % will now get an enclosed shape when plotted

if drawOption==1
    figure('Name','Production Envelope');
    plot(prodEnvVals(:,1),prodEnvVals(:,2),'LineWidth',3)
    xlabel('Growth (/h)')
    ylabel('Target flux (mmol/gDW/h)')
    
end

end

