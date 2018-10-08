function objPFA = objective_PFA(objVARposition,IntIVposition,...
      SRhorizon,SRhorizonIV,B,A,gam,delta)

% objVARpostion is the position of the variable to max SRhorizon impact
% iIntIVposition is the position of the variable we impose penalty for sign
% SRhorizon is the short run horizon of the maximization of the impact of
% objective variable
% B: reduce-form VAR OLS coefficients (nvar*nlags,nvar)
% A: random impact matrix coming from whatever identification, here we use
% one from cholesky is (nvar*nvar).
% gamma: is a (nvar x 1) matrix from the orthogonal matrix

[nvarlags, nvar]         = size(B);
nlags                    = nvarlags/nvar;
IRFs_all                 = zeros(nvar,SRhorizon);
nshocks                  = nvar;

for i_shock = 1:nshocks
      shocks              = zeros(nvar,1);
      shocks(i_shock,1)   = 1;
      % Initialize:
      IRFs_all(:,1,i_shock)   = A*shocks;
      F                       = [IRFs_all(:,1,i_shock)' zeros(1,(nlags-1)*nvar)];
      % Generate IRFs
      for k = 2:SRhorizon
            IRFs_all(:,k,i_shock)    = F*B;
            F                        = [IRFs_all(:,k,i_shock)' F(1:end-nvar)];
      end
end

% Total fluctuations of objVAR
DEN = sum(sum(IRFs_all(objVARposition,:,:).^2)); %independent from gamma

% Initialize:
obj_IRFs(:,1)            = A*gam;
F                        = [obj_IRFs(:,1)' zeros(1,(nlags-1)*nvar)];
% Generate IRFs
for k = 2:SRhorizon
      obj_IRFs(:,k)        = F*B;
      F                    = [obj_IRFs(:,k)' F(1:end-nvar)];
end

if SRhorizonIV == 1
      % Internal instrument
      IVimpact = A*gam; % impact vector
      IVimpact = IVimpact(IntIVposition); % impact on internal instrument
      
      % Objective function
      FEV_objVAR = sum(obj_IRFs(objVARposition,:)); % From IRF to FEV of objVAR
      FEV_objVAR = FEV_objVAR/DEN; % Normalization of FEV
      
      % Penalty Function
      objPFA = FEV_objVAR + delta*IVimpact; % be careful on the sign of delta! ...
      %It defines the sign restriction penalty!
      
      % Negative for fmincon
      objPFA = - objPFA; %BE CAREFUL! IT IS ALREADY NEGATIVE FOR FMINCON!!!
else
      % Internal instrument FEV building
      % Initialize:
      obj_IRFs(:,1)            = A*gam;
      F                        = [obj_IRFs(:,1)' zeros(1,(nlags-1)*nvar)];
      for k = 2:SRhorizonIV
            obj_IRFsIV(:,k)        = F*B;
            F                      = [obj_IRFsIV(:,k)' F(1:end-nvar)];
      end
      FEV_IV = sum(obj_IRFs(IntIVposition,:)); % From IRF to FEV of IV
      FEV_IV = FEV_IV/DEN; % Normalization of FEV
      
      % Objective function
      FEV_objVAR = sum(obj_IRFs(objVARposition,:)); % From IRF to FEV of objVAR
      FEV_objVAR = FEV_objVAR/DEN; % Normalization of FEV
      
      % Penalty Function
      objPFA = FEV_objVAR + delta*FEV_IV; % be careful on the sign of delta! ...
      %It defines the sign restriction penalty!
      
      % Negative for fmincon
      objPFA = - objPFA; %BE CAREFUL! IT IS ALREADY NEGATIVE FOR FMINCON!!!
end


end