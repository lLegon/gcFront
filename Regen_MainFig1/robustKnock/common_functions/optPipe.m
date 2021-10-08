function [screening_results_min, fluxSolution_min] = ...
                optPipe(cobra_model,BM_rxn,target_rxn,max_KOs,branchId, timeLimitRK, endTol, userRxns)
% Modified by Laurence Legon (30/09/21) to only run robustknock, to
% allow termination after a time limit exceeded or after a GC design
% discovered, and to only test a list of user-defined reactions (so essential
% reactions aren't considered when finding KO solutions).
%
%function [optknock_results,robokod_results,optgene_results,screening_results_max,screening_results_min] = optPipe(cobra_model,BM_rxn,target_rxn,max_KOs,branchId)
cwd = cd;
result_dir = [cwd filesep 'results' filesep branchId];
if ~exist(result_dir, 'dir')
    mkdir(result_dir)
end

hasgraphics = usejava('jvm') && usejava('Desktop')&& feature('ShowFigureWindows')

%% Maximum target and biomass production
cobra_model=changeObjective(cobra_model,target_rxn);
FBA.target=optimizeCbModel(cobra_model);

cobra_model=changeObjective(cobra_model,BM_rxn);
FBA.growth=optimizeCbModel(cobra_model);
maxBM=FBA.growth.f;

%clc;
%% Pre-processing
disp('---------Pre-processing----------');

filename = [result_dir filesep branchId '_preprocessing.xls'];

hasrelevant = false;
% % removed because we want to compare initial runs of each algorithm
% if exist(filename)
%     try
%         [~, relevantRxns] = xlsread(filename,'relevantRxns');
%         selectedRxnIDs = find(ismember(cobra_model.rxns, relevantRxns));
%         hasrelevant = true;
%     catch me
%     end
% end

if ~hasrelevant
    [irrel relevantRxns selectedRxnIDs] = irrelevantRxns2(cobra_model,BM_rxn, userRxns);

%     xlwrite(filename, cobra_model.rxns(irrel.lethalRxns) ,'lethalRxns');
%     xlwrite(filename, cobra_model.rxns(irrel.tr_synthRxns) ,'transport-synthRxns');
%     xlwrite(filename, cobra_model.rxns(irrel.blockedRxns) ,'blockedRxns');
%     xlwrite(filename, (irrel.names) ,'all');
%     xlwrite(filename, relevantRxns, 'relevantRxns');
end

% %% OptKnock
% disp('---------Optknock----------');
% 
% [optknock_results]=runOptknock(cobra_model,selectedRxnIDs,BM_rxn,target_rxn,max_KOs,0.5*maxBM);
% 
% if sum(sum(~cellfun(@isempty,optknock_results)))==0
%     disp('No KOs predicted by optknock')
% end
% 
% %TODO: extend optknock result for 3...
% filename = [result_dir filesep branchId '_optknock.xls'];
% xlwrite(filename, optknock_results);
% 
% fprintf('\n')
% %% Robokod
% disp('---------Robokod----------');
% 
% try
% robokod_results=robokod(cobra_model, BM_rxn, target_rxn, max_KOs, selectedRxnIDs);
% catch
%     disp('No KOs predicted by Robokod')
%     robokod_results = cell(1,max_KOs);
% end
% 
% filename = [result_dir filesep branchId '_robokod.xls'];
% xlwrite(filename, robokod_results);
% 
% fprintf('\n')
% %% OptGene
% disp('---------OptGene----------');
% 
% % Change model so that the "irrelevant reactions" don't have a gene association
% % This is needed for OptGene?
% out_model=cobra_model;
% for i=1:length(selectedRxnIDs)
%     in=selectedRxnIDs(i);
%     if in <= length(cobra_model.grRules)
%         out_model.grRules{in}='';
%         out_model.rxnGeneMat(in,:)=0;
%     end
% end
% 
% %Export model
% %TODO
% % Uncomment! this if you would like to export the model, needs sbml binding
% %writeCbModel(cobra_model,'sbml',strcat(result_dir, filesep, branchId, '_model'));
% disp('model created in results folder');
% 
% 
% %Run optgene externally
% disp('Run OptGene in OptFlux for gene deletion:');
% disp('MOMA, 5000 eval, max 3 KOs, SPEA2, Gene Del, maximize BPCY with substrate as coumarate');
% 
% 
% %Read optgene results
% disp('Insert xlsx file with optgene results in the results directory') %'Example: 'rxn1|rxn2|rxn3');
% %Ask for user input( if the results for optgene are already in the folder)
% %to carry on
% if hasgraphics 
%     uiwait(msgbox('Is the optgene_<target>.xlsx file already in the results folder?'));
% else
%     disp('Is the optgene_<target>.xlsx file already in the results folder?');
%     disp('Press any key to continue');
%     pause;
% end
% 
% try
%     [num,optgene_results,raw] = xlsread(strcat(result_dir, filesep, 'optgene_', branchId, '.xlsx'));
% 
%     %if strfind(str,pattern)
%     %Transform genes in reactions
%     optgene_results_rxns= cell(0);
%     for i=1:size(optgene_results,1)
%         [~, tmp]= findRxnsFromGenes(cobra_model, optgene_results(i,:),'listresults',true);
%         if ~isempty(tmp)
%             tmp = tmp(:,1);
%             optgene_results_rxns(end+1,1:length(tmp)) = tmp';
%         end
%     end
%     optgene_results_rxns(find(cellfun(@isempty,optgene_results_rxns))) = {''};
%     optgene_results = optgene_results_rxns;
% catch
%     disp('No KOs predictions by optgene')
%     optgene_results = cell(1,max_KOs);
% end
% 
% %% Enumeration methods
% 
% %OptKnock target: What is the maximum Target when Maximizing Biomass?
% 
% if hasgraphics 
%     uiwait(msgbox('Enumeration can take long time, make sure that you have activated parallel toolbox'));
% else
%     disp('Enumeration can take long time, make sure that you have activated parallel toolbox');
%     disp('Press any key to continue');
%     pause;
% end
% disp('---------Enumeration with OptKnock objective----------');
% 
% filename = [result_dir filesep branchId '_optknock_enumeration.xls'];
% 
% hasoptkonckenum = false;
% 
% if exist(filename)
%     try
%         [~, screening_results_max] = xlsread(filename);
%         %max_screening_results = find(ismember(relevantRxns, cobra_model.rxns));
%         %screening_results_max = cobra_model.rxns(max_screening_results);
%         hasoptkonckenum = true;
%     catch me
%     end
% end
% 
% try
% 
% if ~hasoptkonckenum
%     % Setting the objective for target
%     p = zeros(size(cobra_model.c));
%     targetRxnInd = find(ismember(cobra_model.rxns,  target_rxn));
%     p(targetRxnInd) = -1;
%     evalf = @(model) robustKnockSolution(model,p);
% 
%     % Iterate through all combinations of deletions
%     [max_screening_results,fluxSolution_max] = testRxnDeletion(cobra_model, ... 
%         evalf, max_KOs, selectedRxnIDs);
% 
%     screening_results_max = cobra_model.rxns(max_screening_results);
%     xlwrite(filename, screening_results_max);
% end;
% 
% catch
%     disp('No KOs predictions by optimistic prediction')
%     screening_results_max = cell(1,max_KOs);
% end


%RobustKnock target: What is the minimum Target when Maximizing Biomass?

disp('---------Enumeration with RobustKnock objective----------');

filename = [result_dir filesep branchId '_robustknock_enumeration.xls'];

hasrobustkonckenum = false;
if exist(filename)
    try
        [~, screening_results_min] = xlsread(filename);
        %min_screening_results = find(ismember(relevantRxns, cobra_model.rxns));
        %screening_results_min = cobra_model.rxns(min_screening_results);
        hasrobustkonckenum = true;
    catch me
    end
end

%try

if ~hasrobustkonckenum
    % Setting the objective for target
    p = zeros(size(cobra_model.c));
    targetRxnInd = find(ismember(cobra_model.rxns,  target_rxn));
    p(targetRxnInd) = 1;
    evalf = @(model) robustKnockSolution(model,p);
%changeCobraSolver('ibm_cplex','all',0); % changed to use CPLEX 

    % Iterate through all combinations of deletions
    startTime=tic;
    [min_screening_results,fluxSolution_min] = testRxnDeletion(cobra_model, ... 
        evalf, max_KOs, selectedRxnIDs, timeLimitRK, startTime, endTol);
    screening_results_min = cobra_model.rxns(min_screening_results);
    if toc(startTime)>timeLimitRK
        disp('RobustKnock terminated prematurely')
    end

%    xlwrite(filename, screening_results_min);
end

% catch
%     disp('No KOs predictions by pessimistic prediction')
%     screening_results_min = cell(1,max_KOs);
%     fluxSolution_min=zeros(1,2);
% end

return 
% %{
% %Filter results
% I=find(fluxSolution_min(:,1)>0.5*maxBM);
% screening_i_min=min_screening_results(I,:);
% screening_results_min=cobra_model.rxns(screening_i_min);
% 
% 
% I=find(fluxSolution_max(:,1)>0.5*maxBM);
% screening_i_max=max_screening_results(I,:);
% screening_results_max=cobra_model.rxns(screening_i_max);
% 
% fprintf('\n')
% 
% %}
% %% Organize and analyze results
% 
% % %optknock and optgene results may have only 2 columns
% % if size(optgene_results,2)<max_KOs
% %     optgene_results(:,end+1:max_KOs)={''};
% % end
% % if size(optknock_results,2)<max_KOs
% %     optknock_results(:,end+1:max_KOs)={''};
% % end
% % 
% % %assign method name to each result set
% % optknock_results(:,end+1)={'optknock'};
% % robokod_results(:,end+1)={'robokod'};
% % optgene_results(:,end+1)={'optgene'};
% % screening_results_max(:,end+1)={'screening max'};
% screening_results_min(:,end+1)={'screening min'};
% 
% %join
% %results is for output, all_comb is for use in the next step
% %results=vertcat(optknock_results,robokod_results,optgene_results,screening_results_max,screening_results_min);
% results=vertcat(screening_results_min);
% all_comb=results(:,1:max_KOs);
% 
% %TODO should be sorted!
% %sorting all the rows
% %all_comb = sort(all_comb')';
% 
% %remove empty rows
% a=sum(cellfun(@isempty,all_comb),2);
% all_comb(a==max_KOs,:)=[];
% results(a==max_KOs,:)=[];
% 
% %remove duplicate rows and register duplicates
% tmp = all_comb(:,1);
% for i = 2:max_KOs
%     tmp = strcat(tmp,all_comb(:,i));
% end
% 
% [~,idx]=unique(tmp, 'rows');
% duplicate_i=setdiff(1:length(all_comb),idx);
% duplicates=all_comb(duplicate_i,:);
% all_comb=all_comb(idx,:);
% 
% %remove this part latter
% %all_comb(66,:)=[];
% %all_comb(201,:)=[];
% 
% %Get list of all reactions and analyze (subsystem enrichment)
% [table_subsystems,all_rxns] = analyze_results(all_comb,cobra_model);
% 
% all_rxns(:,2)=printRxnFormula(cobra_model,all_rxns(:,1),false);
% 
% 
% 
% %% Rank combinations
% 
% disp('---------Ranking the KO combinations----------');
% 
% [FVA.min,FVA.max] = fluxVariability(cobra_model);
% 
% for i=1:length(results)
%     comb=results(i,1:max_KOs);
%     comb(~cellfun('isempty',comb));
%     [biomass(i,1) target_min(i,1) target_max(i,1) distance(i,1)] = validateDeletion(cobra_model,...
%         comb,target_rxn,FVA);
% end
% 
% %Rankproduct
% M = [biomass target_min target_max -distance]
% 
% %repeat RankProduct 10 times to acount for randomness of the ranks
% %attributed to parameters with the same values (e.g. multiple double KOs with
% %biomass=0.5)
% all_Ranks=zeros(length(results),4,10);
% for i=1:10
%     temp_rank=rankoptpipe(M);
%     all_Ranks(:,:,i) = [temp_rank.PG,temp_rank.PE,temp_rank.QG,temp_rank.QE];
% end
% %average rank
% avg_Ranks=mean(all_Ranks,3);
% ranks.PG=avg_Ranks(:,1);
% ranks.PE=avg_Ranks(:,2);
% ranks.QG=avg_Ranks(:,3);
% ranks.QE=avg_Ranks(:,4);
% 
% 
% ranks.PG
% 
% 
% 
% %filter out mutants that have max target=0 or biomass<0.1
% remove_ind=find(target_max==0);  
% 
% %This would be too restrictive
% %remove_ind=[remove_ind; find(biomass<0.5*maxBM)];
% remove_ind=[remove_ind; find(biomass<0.1)];
% 
% biomass(remove_ind)=[];
% target_min(remove_ind)=[];
% target_max(remove_ind)=[];
% distance(remove_ind)=[];
% results(remove_ind,:)=[]; 
% 
% ranks.PG(remove_ind,:)=[]; 
% ranks.PE(remove_ind,:)=[]; 
% ranks.QG(remove_ind,:)=[]; 
% ranks.QE(remove_ind,:)=[]; 
% 
% %choose only mutants that have minimal target >0
% choose_ind=find(target_min>0); 
% 
% best_results=results(choose_ind,:);
% best_biomass=biomass(choose_ind);
% best_target_min=target_min(choose_ind);
% best_target_max=target_max(choose_ind);
% best_distance=distance(choose_ind);
% %best_distance=distance(choose_ind);
% 
% %organize
% n=size(results,1);
% table_header(1:max_KOs) = {''};
% table_header = [ table_header {'method','biomass','minimal target','maximal target','distance to Wt','PG','PE','QG','QE'}];
% table_results = table_header;
% table_results(2:n+1,1:max_KOs+1)=results;
% table_results(2:n+1,max_KOs+2)=num2cell(biomass);
% table_results(2:n+1,max_KOs+3)=num2cell(target_min);
% table_results(2:n+1,max_KOs+4)=num2cell(target_max);
% table_results(2:n+1,max_KOs+5)=num2cell(distance);
% 
% table_results(2:n+1,max_KOs+6)=num2cell(ranks.PG);
% table_results(2:n+1,max_KOs+7)=num2cell(ranks.PE);
% table_results(2:n+1,max_KOs+8)=num2cell(ranks.QG);
% table_results(2:n+1,max_KOs+9)=num2cell(ranks.QE);
% 
% n=size(best_results,1);
% table_best_results = table_header;
% table_best_results(2:n+1,1:max_KOs+1)=best_results;
% table_best_results(2:n+1,max_KOs+2)=num2cell(best_biomass);
% table_best_results(2:n+1,max_KOs+3)=num2cell(best_target_min);
% table_best_results(2:n+1,max_KOs+4)=num2cell(best_target_max);
% table_best_results(2:n+1,max_KOs+5)=num2cell(best_distance);
% 
% fprintf('\n')
% %% Best candidates
% 
% disp('---------best candidates----------');
% fprintf('\n')
% 
% [~,I] = sort(best_distance,'ascend');
% a=best_results(I,:);
% disp('Best candidates, taking into account adaptability:')
% %disp(a(1:5,:))
% disp(a(1:end,:))
% 
% fprintf('\n')
% [~,I] = sort(best_target_max,'descend');
% a=best_results(I,:);
% disp('Best candidates, taking into account maximal production:')
% disp(a(1:end,:))
% 
% fprintf('\n')
% [~,I] = sort(best_target_min,'descend');
% a=best_results(I,:);
% disp('Best candidates, taking into account minimal production:')
% disp(a(1:end,:))
% 
% 
% fprintf('\n')
% 
% %% Output results
% disp('---------Exporting to excel----------');
% 
% cd(result_dir);
% 
% 
% %For windows:
% %xlswrite('fisetin_results',all_rxns,'All KOs') 
% %xlswrite('fisetin_results.xls',table_subsystems,'Subsystems enrichment');
% %xlswrite('fisetin_results.xls',table_all_comb,'All KO combinations');
% %xlswrite('fisetin_results.xls',table_best_comb,'Best KO combinations');
% 
% %For mac:
% %uses an external function from file exchange:
% %http://www.mathworks.com/matlabcentral/fileexchange/37560-xlwrite---export-data-to-excel-from-matlab-on-mac-win
% 
% %This one works for linux as well
% %http://www.mathworks.com/matlabcentral/fileexchange/37560-xlwrite---export-data-to-excel-from-matlab-on-mac-win
% 
% 
% 
% xlwrite([branchId '_results.xls'],all_rxns,'All KO reactions'); %put formulas
% xlwrite([branchId '_results.xls'],table_subsystems,'Subsystems enrichment');
% xlwrite([branchId '_results.xls'],table_results,'All KO combinations');
% xlwrite([branchId '_results.xls'],table_best_results,'Best KO combinations');
% 
% disp('Document written');




