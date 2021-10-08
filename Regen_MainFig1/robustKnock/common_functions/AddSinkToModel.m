function [new_model_r] = AddSinkToModel(new_model)
%function [new_model_r] = AddSinkToModel(new_model)
%Adds Malonyl-CoA sink and export reactions
% Deplicated, should be deleted

new_model_r = addReaction( new_model,'malonylCoA_export',{'malcoa[c]','malcoa[e]'},[-1 1]);

new_model_r = addReaction( new_model_r,'sink_malonylCoA',{'malcoa[e]'},[-1]); 

new_model_r = changeRxnBounds( new_model_r, {'sink_malonylCoA'}, [0], 'l');
new_model_r = changeRxnBounds( new_model_r, {'sink_malonylCoA'}, [1000], 'u');

end

