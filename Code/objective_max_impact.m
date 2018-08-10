function obj_FEV = objective_max_impact(which_variable,A,gam)

% which_variable: It select the variable that the shock is maximizing (TFP)
% A: random impact matrix coming from whatever identification, here we use
% one from cholesky is (nvar*nvar). 
% gamma: is a (nvar x 1) matrix 

obj = A*gam;
obj_FEV = - obj(which_variable);
if obj_FEV > 0
    obj_FEV = - obj_FEV;
end