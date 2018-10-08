function [impact, gam] = identification_GPFA(A,B,SRhorizon,SRhorizonIV,EBPposition,...
      Uposition,IVposition,delta)

% A chol impact matrix (nvar,nvar)
% Reduce-form regressors coefficient (1+nvar*nlags,nvar)
% SRhorizon is which horizon to max news shock
% EBPposition is position of excess bond premium in the VAR
% Uposition is position of JLN macro uncertainty in the VAR
% CFposition is position of normalized cashflow in the VAR
% delta is penalty function parameter (scalar) - same for both

% Technical values
B                     = B(2:end,:); % remove the constant
[nvarlags, nvar]      = size(B);
nlags                 = nvarlags/nvar;

% Defining initial values
D              = eye(nvar);
gamF_zero   = D(:,1); %financial shock impact vector (initial value)
gamU_zero   = D(:,2); %uncertainty shock impact vector (initial value)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Step - Identifying gamF - Financial Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the objective function of Step1
objF = @(gamF) objective_PFA(EBPposition,IVposition,SRhorizon,...
      SRhorizonIV,B,A,gamF,delta); 

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

% %Constraint that News shocks have no contemporaneous effect on TFP
% Me3       = 1; %  Me = no. of equality constraints
% Beq3      = zeros(Me3,1); % Beq is (Me x 1) where
% Aeq3      = zeros(Me3,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
% Aeq3(1,1) = 1; %zero-impact of news on TFP

% Optimization
[gamF_opt] = fmincon(objF, gamF_zero,[],[],[],[],[],[],...
      @(gamF) constraint_orthogonality(gamF),options);
% [gamF_opt] = fmincon(objF, gamF_zero,[],[],Aeq3,Beq3,[],[],...
%       @(gamF) constraint_barskysims(gamF),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (gamF_opt'*gamF_opt - 1)^2 > 10^(-10) %|| gamF_opt(1)^2 > 10^(-12)
      warning('The problem is not consistent with the constraints.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Step - Identifying gamU - Uncertainty Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the objective function of Step2
delta2 = - delta; % Be careful on the sign of delta!
objU = @(gamU) objective_PFA(Uposition,IVposition,SRhorizon,SRhorizonIV,B,A,gamU,delta2);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

%Optimization - Notice the contraint for gamNews and gamTFP
gamU_opt = fmincon(objU, gamU_zero,[],[],[],[],[],[],...
      @(gamU) constraint_orthogonality([gamU]),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (gamU_opt'*gamU_opt - 1)^2 > 10^(-10) || sum(sum(([gamF_opt gamU_opt]'*[gamF_opt gamU_opt] - eye(2)).^2)) > 10^(-10)
      warning('The problem is not consistent with the constraints.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Step - Obtain impact matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gam          = [gamF_opt gamU_opt];
impact1      = A*gam;
impact       = [impact1 zeros(size(impact1,1),size(impact1,1)-size(impact1,2))];

end
