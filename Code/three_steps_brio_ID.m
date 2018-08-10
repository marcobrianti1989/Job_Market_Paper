function [impact, gam] = three_steps_brio_ID(A,B,horizon,TFPposition,VXOposition)

% A chol impact matrix
% Reduce-form regressors coefficient
% horizon is which horizon to max news shock
% TFPposition is position of TFP in the VAR
% VXOposition is position of implied volatility variable in the VAR

% Technical values
B                     = B(2:end,:);
[nvarlags, nvar]      = size(B);
nlags                 = nvarlags/nvar;

% Defining initial values
D              = eye(nvar);
gamNews_zero   = D(:,1); %news shock impact vector (initial value)
gamTFP_zero    = D(:,2); %surprize TFP shock impact vector (initial value)
gamUnc_zero    = D(:,3); %uncertainty shock impact vector (initial value)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Step - Identifying gamNews for news shock - Similar to B&S(2012)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the objective function of Step1
objNews = @(gamNews) objective_barskysims(TFPposition,horizon,B,A,gamNews); % usual is BS

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

%Constraint that News shocks have no contemporaneous effect on TFP
Me3       = 1; %  Me = no. of equality constraints
Beq3      = zeros(Me3,1); % Beq is (Me x 1) where
Aeq3      = zeros(Me3,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
Aeq3(1,1) = 1; %zero-impact of news on TFP

% Optimization
[gamNews_opt] = fmincon(objNews, gamNews_zero,[],[],Aeq3,Beq3,[],[],...
      @(gamNews) constraint_barskysims(gamNews),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (gamNews_opt'*gamNews_opt - 1)^2 > 10^(-10) || gamNews_opt(1)^2 > 10^(-12)
      warning('The problem is not consistent with the constraints.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Step - Identifying gamTFP - Surprise TFP shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the objective function of Step2
objTFP = @(gamTFP) objective_max_impact(TFPposition,A,gamTFP);

%Optimization - Notice the contraint for gamNews and gamTFP
gamTFP_opt = fmincon(objTFP, gamTFP_zero,[],[],[],[],[],[],...
      @(gamTFP) constraint_orthogonality([gamNews_opt gamTFP]),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (gamTFP_opt'*gamTFP_opt - 1)^2 > 10^(-10) || sum(sum(([gamNews_opt gamTFP_opt]'*[gamNews_opt gamTFP_opt] - eye(2)).^2)) > 10^(-10)
      warning('The problem is not consistent with the constraints.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third Step - Identifying gamUnc - Uncertainty shocks - Similar to Basu&Bundick (2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the objective function of Step3
objUnc = @(gamUnc) objective_max_impact(VXOposition,A,gamUnc);

% Optimization
gamUnc_opt = fmincon(objUnc, gamUnc_zero,[],[],[],[],[],[],...
      @(gamUnc) constraint_orthogonality([gamNews_opt gamTFP_opt gamUnc]),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (gamUnc_opt'*gamUnc_opt - 1)^2 > 10^(-10) || sum(sum(([gamNews_opt gamTFP_opt gamUnc_opt]'*[gamNews_opt gamTFP_opt gamUnc_opt] - eye(3)).^2)) > 10^(-10)
      warning('The problem is not consistent with the constraints.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Step - Obtain impact matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gam          = [gamNews_opt gamTFP_opt gamUnc_opt];
impact1      = A*gam;
impact       = [impact1 zeros(size(impact1,1),size(impact1,1)-size(impact1,2))];

end
