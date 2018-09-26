function [impact, gam] = just_news_ID(A,B,horizon,position)

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
objNews = @(gamNews) objective_barskysims(position,horizon,B,A,gamNews); % usual is BS

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
% Last Step - Obtain impact matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gam          = gamNews_opt;
impact1      = A*gam;
impact       = [impact1 zeros(size(impact1,1),size(impact1,1)-size(impact1,2))];

end
