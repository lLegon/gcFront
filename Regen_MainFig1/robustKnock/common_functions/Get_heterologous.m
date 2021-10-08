function [new_model sink] = Get_heterologous(model,pathway)
%function [new_model sink] = Get_heterologous(model,pathway)
% Adds the heterologous pathways for biophenolics
%      INPUT:
%               model: COBRA model to add the pathways
%  	            pathway: one of the four branches 'resveratrol', 'pelargonidin', 'fisetin', 'quercetin'

%modified version 15/12/2015: 
%	*use coumarate as a substrate: don't insert the common branch
%	*in the resveratrol pathway: go only up  to the trans-resveratrol point (remove vv-ROMT reactions)
%	*use fisetin and quercetin pathways from Stalhut, 2015

%NOTES
% IN the fisetin and quercetin pathway there migth more reactions ooccuring in vivo. 
%F3bH, FLS and FMO can provide  an alternative pathway for the production of quercertin, through eriodictyol and taxifolin.
%However, since these alternative pathways are stoichiometrically equivalent, they were ommited for simplicity.

%% Assigning the sink reaction
%If the resveratrol branch is added: optimize 'EX_transresveratrol'
%If the pelargonidin branch is added: optimize 'EX_pelargonidin'
%If the fisetin branch is added: optimize 'EX_fisetin'
%If the quecetin branch is added: optimize 'EX_quercetin'
%Throws error if different
switch pathway
    case 'resveratrol' 
        sink = 'sink_transresveratrol';
    case 'pelargonidin' 
        sink = 'sink_pelargonidin';
    case 'fisetin' 
        sink = 'sink_fisetin';
    case 'quercetin'
        sink = 'sink_quercetin'
    otherwise
        error('Get_heterologous:pathway',...
        ['Pathway for ''' pathway ''' does not exist.\n Pathway should be one of the four branches ''resveratrol'', ''pelargonidin'', ''fisetin'', or ''quercetin''.'])
end

new_model=model;

%% Common branch
% 4-coumarate coA ligase (4CL):4-coumarate + ATP + coA ? 4-coumaryl-coA + AMP + diphosphate
new_model=addReaction(new_model,'4CL',{'4_coumarate[c]','atp[c]','coa[c]','4_coumaryl_coa[c]','amp[c]','ppi[c]'},[-1 -1 -1 1 1 1]);
new_model.rxnNames{end,1}='4-coumarate coA ligase';
new_model.grRules{end,1}='4CL';

%add sinks for coumarate and coumaryl-coA
new_model=addReaction(new_model,'sink_coumarate',{'4_coumarate[c]'},[-1]);
new_model.rxnNames{end,1}='sink coumarate';
new_model.lb(end)=0;
new_model=addReaction(new_model,'sink_coumarylcoA',{'4_coumaryl_coa[c]'},[-1]);
new_model.rxnNames{end,1}='sink coumaryl-coA';
new_model.lb(end)=0;

disp('Common reaction 4CL added');
new_model.description='extended with the common branch of the polyphenolics'




%% RESVERATROL branch
if strcmp(pathway,'resveratrol')==1
%Stilbene synthase(STS):4-coumaryl-coA + 3 malonyl-coA + 3H+ ? trans-resveratrol (stilbene) + 4 CO2 + 4 coA
    new_model=addReaction(new_model,'STS',{'4_coumaryl_coa[c]','malcoa[c]','h[c]','trans_resveratrol[c]','co2[c]','coa[c]'},[-1 -3 -3 1 4 4]);
    new_model.rxnNames{end,1}='stilbene synthase';
    new_model.grRules{end,1}='STS';

%add export and sink for trans-resveratrol 
    new_model=addReaction(new_model,'EX_transresveratrol',{'trans_resveratrol[c]','trans_resveratrol[e]'},[-1 1]);
    new_model.rxnNames{end,1}='export trans-resveratrol';
    new_model=addReaction(new_model,'sink_transresveratrol',{'trans_resveratrol[e]'},[-1]);
    new_model.rxnNames{end,1}='sink transresveratrol';
    new_model.lb(end)=0;
    
    disp('Resveratrol branch added');
    new_model.description='extended with resveratrol branch';
end
    
    
%% Pelargonidin branch
if strcmp(pathway,'pelargonidin')==1

% Naringenin- chalcone synthase(CHS): 4-coumaryl-coA + 3 malonyl-coA + 3 H+ ? naringenin chalcone + 3 CO2 + 4 coA
    new_model=addReaction(new_model,'CHS',{'4_coumaryl_coa[c]','malcoa[c]','h[c]','naringenin_chalcone[c]','co2[c]','coa[c]'},[-1 -3 -3 1 3 4]);
    new_model.rxnNames{end,1}='naringenin_chalcone synthase';
    new_model.grRules{end,1}='CHS';

% Chalcone isomerase(CHI): naringenin chalcone ? naringenin
    new_model=addReaction(new_model,'CHI',{'naringenin_chalcone[c]','naringenin[c]'},[-1 1]);
    new_model.rxnNames{end,1}='chalcone isomerase';
    new_model.grRules{end,1}='CHI';

% Flavonone 3-dioxygenase(F3H): Naringenin + 2-oxoglutarate + O2 ? dihydrokaempferol + succinate + CO2  IRREVERSIBLE
    new_model=addReaction(new_model,'F3H',{'naringenin[c]','akg[c]','o2[c]','dihydrokaempferol[c]','suc[c]','co2[c]'},[-1 -1 -1 1 1 1],'false');
    new_model.rxnNames{end,1}='flavonone 3-dioxygenase';
    new_model.grRules{end,1}='F3H';

% Dihydroflavonol 4-reductase(DFR): Dihydroxykaempferol + NAPH + H+ ? Leucopelargonidin + NAP+ 
    new_model=addReaction(new_model,'DFR',{'dihydrokaempferol[c]','nadph[c]','h[c]','leucopelargonidin[c]','nadp[c]'},[-1 -1 -1 1 1]);
    new_model.rxnNames{end,1}='dihydroflavonol 4-reductase';
    new_model.grRules{end,1}='DFR';

% Anthocyanidin synthase(ANS): Leucopelargonidin + 2-oxoglutarate + O2 ? pelargonidin + succinate + CO2 + H+ + 2 H2O IRREVERSIBLE
    new_model=addReaction(new_model,'ANS',{'leucopelargonidin[c]','akg[c]','o2[c]','pelargonidin[c]','suc[c]','co2[c]','h[c]','h2o[c]'},[-1 -1 -1 1 1 1 1 2],'false');
    new_model.rxnNames{end,1}='anthocyanidin synthase';
    new_model.grRules{end,1}= 'ANS';

%add sinks for naringenin_chalcone, naringenin, dihydrokaempferol and
%leucopelargonidin
    new_model=addReaction(new_model,'sink_naringenin_chalcone',{'naringenin_chalcone[c]'},[-1]);
    new_model.rxnNames{end,1}='sink naringenin chalcone';
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_naringenin',{'naringenin[c]'},[-1]);
    new_model.rxnNames{end,1}='sink naringenin';
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_dihydrokaempferol',{'dihydrokaempferol[c]'},[-1]);
    new_model.rxnNames{end,1}='sink dihydrokaempferol';
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_leucopelargonidin',{'leucopelargonidin[c]'},[-1]);
    new_model.rxnNames{end,1}='sink leucopelargonidin';
    new_model.lb(end)=0;

%add export and sink for pelargonidin (target compound)
    new_model=addReaction(new_model,'EX_pelargonidin',{'pelargonidin[c]','pelargonidin[e]'},[-1 1]);
    new_model.rxnNames{end,1}='export pelargonidin';
    new_model=addReaction(new_model,'sink_pelargonidin',{'pelargonidin[e]'},[-1]);
    new_model.rxnNames{end,1}='sink pelargonidin';
    new_model.lb(end)=0;
    
    disp('Pelargonidin branch added');
    new_model.description='extended with pelargonidin branch';
end


%% Quercetin branch
if strcmp(pathway,'quercetin')==1

% Naringenin- chalcone synthase(CHS): 4-coumaryl-coA + 3 malonyl-coA + 3 H+ ? naringenin chalcone + 3 CO2 + 4 coA
    new_model=addReaction(new_model,'CHS',{'4_coumaryl_coa[c]','malcoa[c]','h[c]','naringenin_chalcone[c]','co2[c]','coa[c]'},[-1 -3 -3 1 3 4]);
    new_model.rxnNames{end,1}='naringenin_chalcone synthase';
    new_model.grRules{end,1}='CHS';

% Chalcone isomerase(CHI): naringenin chalcone ? naringenin
    new_model=addReaction(new_model,'CHI',{'naringenin_chalcone[c]','naringenin[c]'},[-1 1]);
    new_model.rxnNames{end,1}='chalcone isomerase';
    new_model.grRules{end,1}='CHI';

%Flavonoid 3-hydroxylase (F3H): Naringenin + O2 +alfa-oxoglutarate --> Dihydrokaempferol + CO2 + succinate
    new_model=addReaction(new_model,'F3H',{'naringenin[c]','akg[c]','o2[c]','dihydrokaempferol[c]','suc[c]','co2[c]'},[-1 -1 -1 1 1 1]);
    new_model.rxnNames{end,1}='flavonone 3-hydroxylase';
    new_model.grRules{end,1}='F3H';

%Flavonol synthase (FLS) (reversible): Dihydrokaempferol + 2-oxoglutarate +O2 --> kaempferol + succinate + CO2  IRREVERSIBLE
    new_model=addReaction(new_model,'FLS',{'dihydrokaempferol[c]','akg[c]','o2[c]','kaempferol[c]','suc[c]','co2[c]'},[-1 -1 -1 1 1 1],'false');
    new_model.rxnNames{end,1}='Flavonol synthase';
    new_model.grRules{end,1}='FLS';

%Flavonoid 3'Monooxygenase-reductase(FMO-R):kaempferol + nadph + O2 ? quercetin + nadp + h2o IRREVERSIBLE
    new_model=addReaction(new_model,'FMO',{'kaempferol[c]','nadph[c]','o2[c]','quercetin[c]','nadp[c]','h2o[c]'},[-1 -1 -1 1 1 1],'false');
    new_model.rxnNames{end,1}='Flavonoid 3''monooxygenase';
    new_model.grRules{end,1}='FMO';

%add sinks for naringenin chalcone, naringenin, dihydrokaempferol and kaempferol
    new_model=addReaction(new_model,'sink_naringenin_chalcone',{'naringenin_chalcone[c]'},[-1]);
    new_model.rxnNames{end,1}='sink naringenin chalcone';   
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_naringenin',{'naringenin[c]'},[-1]);
    new_model.rxnNames{end,1}='sink naringenin';   
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_dihydrokaempferol',{'dihydrokaempferol[c]'},[-1]);
    new_model.rxnNames{end,1}='sink dihydrokaempferol';   
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_kaempferol',{'kaempferol[c]'},[-1]);
    new_model.rxnNames{end,1}='sink kaempferol';   
    new_model.lb(end)=0;   

%add export and sink for quercetin (target compound)
    new_model=addReaction(new_model,'EX_quercetin',{'quercetin[c]','quercetin[e]'},[-1 1]);
    new_model.rxnNames{end,1}='export quercetin '; 
    new_model=addReaction(new_model,'sink_quercetin',{'quercetin[e]'},[-1]);
    new_model.rxnNames{end,1}='sink quercetin'; 
    new_model.lb(end)=0;

    disp('Quercetin branch added');
    new_model.description='extended with fisetin branch';
end


%% Fisetin branch
%closely related to the quercetin branch
if strcmp(pathway,'fisetin')==1
%Naringenin- chalcone synthase + chalcone reductase(CHS-R): 4-coumaryl-coA + 3 malonyl-coA + 3 H+ ? isoliquiritigenin + 3 CO2 + 4 coA
    new_model=addReaction(new_model,'CHSR',{'4_coumaryl_coa[c]','malcoa[c]','nadph[c]','isoliquiritigenin[c]','coa[c]','nadp[c]'},[-1 -3 -1 1 4 1]);
    new_model.rxnNames{end,1}='chalcone synthase reductase';
    new_model.grRules{end,1}= 'CHS_R';

% Chalcone isomerase(CHI): isoliquiritigenin --> liquiritigenin
    new_model=addReaction(new_model,'CHI',{'isoliquiritigenin[c]','liquiritigenin[c]'},[-1 1]);
    new_model.rxnNames{end,1}='chalcone isomerase';
    new_model.grRules{end,1}='CHI';

%Flavonoid 3-hydroxylase (F3H): Liquiritigenin + O2 +alfa-oxoglutarate --> Garbanzol + CO2 + succinate
    new_model=addReaction(new_model,'F3H',{'liquiritigenin[c]','akg[c]','o2[c]','garbanzol[c]','suc[c]','co2[c]'},[-1 -1 -1 1 1 1]);
    new_model.rxnNames{end,1}='flavonone 3-hydroxylase';
    new_model.grRules{end,1}='F3H';

%Flavonol synthase (FLS) (reversible): Garbanzol + 2-oxoglutarate +O2 --> resokaempferol + succinate + CO2  IRREVERSIBLE
    new_model=addReaction(new_model,'FLS',{'garbanzol[c]','akg[c]','o2[c]','resokaempferol[c]','suc[c]','co2[c]'},[-1 -1 -1 1 1 1],'false');
    new_model.rxnNames{end,1}='Flavonol synthase';
    new_model.grRules{end,1}='FLS';
    
%Flavonoid 3'Monooxygenase reductase(FMO-R):Resokaempferol + nadph + O2 ? fisetin + nadp + h2o IRREVERSIBLE
new_model=addReaction(new_model,'FMO',{'resokaempferol[c]','nadph[c]','o2[c]','fisetin[c]','nadp[c]','h2o[c]'},[-1 -1 -1 1 1 1],'false');
new_model.rxnNames{end,1}='Flavonoid 3''monooxygenase';
new_model.grRules{end,1}='FMO';

%add sinks for isoliquiritigenin, liquiritigenin, garbanzol and resokaempferol
    new_model=addReaction(new_model,'sink_isoliquiritigenin',{'isoliquiritigenin[c]'},[-1]);
    new_model.rxnNames{end,1}='sink isoliquiritigenin';   
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_liquiritigenin',{'liquiritigenin[c]'},[-1]);
    new_model.rxnNames{end,1}='sink liquiritigenin';   
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_garbanzol',{'garbanzol[c]'},[-1]);
    new_model.rxnNames{end,1}='sink garbanzol';   
    new_model.lb(end)=0;
    new_model=addReaction(new_model,'sink_resokaempferol',{'resokaempferol[c]'},[-1]);
    new_model.rxnNames{end,1}='sink resokaempferol';   
    new_model.lb(end)=0;

%add export and sink for fisetin (target compound)
    new_model=addReaction(new_model,'EX_fisetin',{'fisetin[c]','fisetin[e]'},[-1 1]);
    new_model.rxnNames{end,1}='export fisetin '; 
    new_model=addReaction(new_model,'sink_fisetin',{'fisetin[e]'},[-1]);
    new_model.rxnNames{end,1}='sink fisetin'; 
    new_model.lb(end)=0;
    
    disp('Fisetin branch added');
    new_model.description='extended with fisetin branch';
 
end



%% Add subsystem fields
for i=length(model.rxns)+1:length(new_model.rxns)
    new_model.subSystems{i,1}='heterologous';
end


end

