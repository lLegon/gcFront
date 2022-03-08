function state=myPlotPareto3obj(defaultFitness,wtEnv,fitnessLimit,options,state,~,~)
% plot a pareto front of the population during multiobjective GA

%tic

scores=sortrows(-state.Score,3);

if max(scores(:,3))<=0
    % not coupled- plot coupling strength
    
    if state.Generation==0
        % label plot
        %title('No coupled designs found yet')
        xlabel('Growth rate (/h)')
        ylabel('Coupling strength')
        xlim([0,inf])
        ylim([min(setdiff(scores(:,3),defaultFitness)),0])    
    end
    
    scatter(scores(:,1),scores(:,3),50,'kx','LineWidth',3)
    
else
    % change product values of uncoupled designs in scores so they will not appear
    % on the plot of the pareto front
    scores(scores(:,3)<=0,2)=-1;
    
    % get the scatter plot points
    currPlot=findobj(get(gca,'Children'));
    
    % colour points by their coupling strength
    colPoints=[0,0.01,1;
        0.55,0.8,0;
        0.3,0.9,0.6];
    
    colPoints2=[colPoints(1,1):(colPoints(1,2)-colPoints(1,1))/30:colPoints(1,2), colPoints(1,2):(colPoints(1,3)-colPoints(1,2))/30:colPoints(1,3);
        colPoints(2,1):(colPoints(2,2)-colPoints(2,1))/30:colPoints(2,2), colPoints(2,2):(colPoints(2,3)-colPoints(2,2))/30:colPoints(2,3);
        colPoints(3,1):(colPoints(3,2)-colPoints(3,1))/30:colPoints(3,2), colPoints(3,2):(colPoints(3,3)-colPoints(3,2))/30:colPoints(3,3);]';
    
    cols=colormap(colPoints2);
    
    colInds=ceil( scores(:,3) * size(cols,1)/2 );
    uncoupled=colInds<=0;
    colInds(uncoupled)=1;
    
    currCols=cols(colInds,:);
    currCols(uncoupled,:)=1;
    
    if length(currPlot)==2
        % the production envelope and a previous generation of scatter
        % points have been plotted- can just adjust position and color of
        % existing points

        
        currPlot(1).XData=scores(:,1);
        currPlot(1).YData=scores(:,2);
        currPlot(1).CData=currCols;
        
        
    else
        % coupling has just been discovered- plot production envelope and
        % show designs against production envelope
        
        disp(['Coupled strain discovered- GA has been running for ', num2str(toc(state.StartTime)), ' seconds'])
        
        %title('Coupled')
        xlabel('Growth rate (/h)')
        ylabel('Target flux (mmol/gDW/h)')
        

        % removing growth values for uncoupled designs in state.Score-
        % otherwise uncoupled designs with high growth rates will persist
        state.Score(state.Score(:,3)>=0,1)=0;
        
        plot(wtEnv(:,1),wtEnv(:,2),'k','LineWidth',3)
        hold on
        scatter(scores(:,1),scores(:,2),50,currCols,'x','LineWidth',3)
        hold off
        
        xlim([0,ceil(max(wtEnv(:,1)))])
        ylim([0,ceil(max(wtEnv(:,2)))])
        colbar=colorbar;
        colbar.TickLabelsMode='manual';
        colbar.Ticks=0:0.25:1;
        colbar.TickLabels=num2str( (0:0.5:2)' );
        colbar.Label.String='Coupling strength';
        
    end

    % check if fitness limit exceeded
    if max(scores(:,2))>fitnessLimit
        disp('Algorithm terminating- fitness limit exceeded')
        state.StopFlag=1;
    end
    
end

% find % of time/generations that has expired
gensUsed=100*state.Generation/options.MaxGenerations;
timeUsed=100*toc(state.StartTime) / options.MaxTime;
% round to 1 dp and display as title of plot
if gensUsed>timeUsed
    gensUsed=round(gensUsed*10)/10;
    title([num2str(gensUsed), '% of maximum generations completed'])
elseif timeUsed>gensUsed
    timeUsed=round(timeUsed*10)/10;
    title([num2str(timeUsed), '% of maximum time used'])
else
    timeUsed=round(timeUsed*10)/10;
    title([num2str(timeUsed), '% of maximum generations/maximum time used'])
end


% toc
% disp(' ')


end