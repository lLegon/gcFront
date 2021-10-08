function [ Constrained_model ] = Insert_constraints( Model,orig_constraints)
%function [ Constrained_model ] = Insert_constraints( Model,orig_constraints)
% Input: Model: SBML model, 
%        orig_constraints: Include the original constraints if true, and +-20% tolerance if false or non-existent
% Output:Constrained_model: SBML model including constraints

%% constraints for wild type provided by Elisabeth Zelle 
if(exist('orig_constraints')&&orig_constraints)

    Model = changeRxnBounds( Model ,'EX_glc(e)',-4.98,'u');
    Model = changeRxnBounds( Model ,'EX_glc(e)',-5.796,'l');

    Model = changeRxnBounds( Model ,'EX_o2(e)',-3.9,'u');
    Model = changeRxnBounds( Model ,'EX_o2(e)',-4.38,'l');

    Model = changeRxnBounds( Model ,'AC_d',-(0.69/15),'u');
    Model = changeRxnBounds( Model ,'AC_d',-(0.79/15),'l');

    Model = changeRxnBounds( Model ,'SUC_d',(1.01/15),'u');
    Model = changeRxnBounds( Model ,'SUC_d',(0.47/15),'l');

    Model = changeRxnBounds( Model ,'EX_lac_L(e)',(18/25),'u');%
    Model = changeRxnBounds( Model ,'EX_lac_L(e)',(8/25),'l');%

    Constrained_model = changeRxnBounds( Model ,'ACC_glycogen_c',-15,'l');%

    %Constrained_model = changeRxnBounds( Constrained_model, {'EX_glc(e)'}, -5.4, 'b');
    return;
end

    Model = changeRxnBounds( Model ,'EX_glc(e)',-4.98*0.8,'u');%-3.984
    Model = changeRxnBounds( Model ,'EX_glc(e)',-5.796*1.2,'l');%-6.9552

    Model = changeRxnBounds( Model ,'EX_o2(e)',-3.9*0.8,'u');%-3.12
    Model = changeRxnBounds( Model ,'EX_o2(e)',-4.38*1.2,'l');%-5.26

    Model = changeRxnBounds( Model ,'AC_d',-(0.69/15)*0.8,'u');%-0.0368
    Model = changeRxnBounds( Model ,'AC_d',-(0.79/15)*1.2,'l');%-0.632

    Model = changeRxnBounds( Model ,'SUC_d',(1.01/15)*1.2,'u');%0.808
    Model = changeRxnBounds( Model ,'SUC_d',(0.47/15)*0.8,'l');%0.2506667

    Model = changeRxnBounds( Model ,'EX_lac_L(e)',(18/25)*1.2,'u');%0.864
    Model = changeRxnBounds( Model ,'EX_lac_L(e)',(8/25)*0.8,'l');%0.256

    Constrained_model = changeRxnBounds( Model ,'ACC_glycogen_c',-15*1.2,'l');%18
    return;
end

