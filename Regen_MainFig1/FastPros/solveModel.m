function [sol, objv, status, lambda]= solveModel(model)
%%solveModel  モデルを解く（モデルに対してFBAを実行）

[sol, objv, status, extra] = glpk(model.c,model.S,model.b,model.lb,model.ub,...
        repmat('S',size(model.S,1),1),repmat('C',size(model.S,2),1),-1);
lambda = extra.lambda;