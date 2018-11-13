clear
close all

% x  = rho*x(-1) + rh0*y(-1) + rho*z(-2) + epsU + epsF;
% y  = rho*x(-1) + rh0*y(-1) + rho*z(-2) + epsU + epsF;
% z  = rho*x(-1) + rh0*y(-1) + rho*z(-2) + epsU - epsF;

% Parameterization
rho1    = 0.8;
rho2    = 0.1;
rho3    = 0.1;
rho4    = 0.9;
sigU    = 1;
sigF    = 1;
alp     = 5;
A11     = 1;
A21     = 1;
A12     = 0.1;%alp*A11;
A22     = 0.1; %alp*A21;
% Number of Observations
T       = 100000;

% Predetermined matrices
U       = zeros(T,1);
F       = zeros(T,1);
C       = zeros(T,1);
shocksU = sigU*randn(T,1);
shocksF = sigF*randn(T,1);

% Structural Economy
for i = 2:T 
      U(i) = rho1*U(i-1) + rho2*C(i-1) + A11*shocksU(i) + A12*shocksF(i);
      C(i) = rho3*U(i-1) + rho4*C(i-1) + A21*shocksU(i) - A22*shocksF(i);   
end

% Define the system1
system_names  = {'U','C'};
for i = 1:length(system_names)  
      system(:,i) = eval(system_names{i});   
end
Uposition   = find(strcmp('U', system_names));
IVposition  = find(strcmp('C', system_names));

% Cholesky decomposition
nlags           = 1;
[A,B,res,sigma] = sr_var(system, nlags);
corr_res        = corr(res);


fprintf(1, '\n');  
fprintf(1, '\n');  
fprintf(1, '\n');  
disp(sprintf('Covariance of residuals is          : %f',sigma(1,2)));
fprintf(1, '\n'); 
disp(sprintf('Correlation matrix of residuals is  : %f',corr_res(1,2)));
% Generalized Penalty Function Approach
objgam     = 1;
delta      = 0;
threshold  = -1;
add        = 0.01;
j = 1;
while objgam >= -1 %10^(-4)
      while objgam >= threshold
            warning off
            SRhorizon       = 1;
            SRhorizonIV     = 1;
            [impact, gamma] = identification_GPFA(A,B,SRhorizon,...
                  SRhorizonIV,Uposition,...
                  Uposition,IVposition,-delta);
            gamF   = gamma(:,1);
            gamU   = gamma(:,2);
            gamU(1);
            gamU1(j)  = gamU(1);
            gamU2(j)  = gamU(2);
            gamF1(j)  = gamF(1);
            gamF2(j)  = gamF(2);
            gamF'*gamU
            objgam(j) = gamF'*gamU;
            delta  = delta + add;
            j = j + 1;
            %pause(0.005)
      end
      threshold = threshold/10;
      add       = add/10;
end
figure
subplot(2,2,1)
plot(gamU1)
subplot(2,2,2)
plot(gamF1)
subplot(2,2,3)
plot(gamU2)
subplot(2,2,4)
plot(gamF2)

figure
plot(objgam)


fprintf(1, '\n');  
fprintf(1, '\n');  
fprintf(1, '\n');  
disp(sprintf('Correlation between shocks is: %f',objgam));
fprintf(1, '\n');  
delta = - delta; % Negative because the first variable is EBP and cash respons negatively on impact

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 5;
which_correction  = 'none';
blocksize         = 4;
[beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian(B,nburn,res,nsimul,which_correction,blocksize);

% Get "bootstrapped A" nsimul times
for i_simul = 1:nsimul    
      % Cholesky decomposition     
      [A_boot(:,:,i_simul),B_boot(:,:,i_simul),~,~] = ...
            sr_var(data_boot2(:,:,i_simul), nlags);     
      % GPFA identification strategy    
      warning off    
      [impact_boot(:,:,i_simul), gamma_boot(:,:,i_simul)] = ...
            identification_GPFA(A_boot(:,:,i_simul),B_boot(:,:,i_simul),...
            SRhorizon,SRhorizonIV,Uposition,Uposition,IVposition,delta);   
      i_simul;   
end

% Generate IRFs with upper and lower bounds
sig1                       = 0.05;
sig2                       = 0.025;
H                          = 20;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(impact,impact_boot,B,B_boot,H,sig1,sig2);

% Theoretical IRF - U shock
xtheory     = zeros(H,2);
ytheory      = zeros(H,2);
xtheory(1,2) = A11*sigU;
ytheory(1,2) = A21*sigU;
for i = 2:H    
      xtheory(i,2) = rho1*xtheory(i-1,2) + rho2*ytheory(i-1,2);
      ytheory(i,2) = rho3*xtheory(i-1,2) + rho4*ytheory(i-1,2);      
end

% Theoretical IRF - F shock
xtheory(1,1) =   A12*sigF;
ytheory(1,1) = - A22*sigF;
for i = 2:H  
      xtheory(i,1) = rho1*xtheory(i-1,1) + rho2*ytheory(i-1,1);
      ytheory(i,1) = rho3*xtheory(i-1,1) + rho4*ytheory(i-1,1);
end
theoryIRF = zeros(2,H,2);
theoryIRF(1,:,:) = xtheory;
theoryIRF(2,:,:) = ytheory;

% Create and Printing figures for IRFs
base_path         = pwd;
which_ID          = 'GPFA_Compustat_3lags_gamgamZero_GDPDef';
use_current_time  = 1; % (don't) save the time
which_shocks      = [1 2]; %[Uposition];
shocknames        = {'Financial Shock','Uncertainty Shock'};
print_figs = 'no';
plot_IRFs_Empirical_Theoretical_2CIs(IRFs,ub1,lb1,ub2,lb2,theoryIRF,H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)

% plot_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
%       system_names,which_ID,print_figs,use_current_time,base_path)

% Get Structural Shocks
ss = get_structural_shocks_general(A,gamma,res,which_shocks);
ssF = ss(:,1);
ssU = ss(:,2);
corrF = corr(ssF,shocksF(2:end));
corrU = corr(ssU,shocksU(2:end));
fprintf(1, '\n');  
fprintf(1, '\n');  
fprintf(1, '\n');  
disp(sprintf('Correlation between true U shocks and estimated ones is: %f',corrU));
fprintf(1, '\n');  
disp(sprintf('Correlation between true F shocks and estimated ones is: %f',corrF));
fprintf(1, '\n');  
fprintf(1, '\n'); 
fprintf(1, '\n'); 
disp('Stop code here.')
return

% Get variance Decomposition
N = null(gamma');
D_null = [gamma N];
impact_vardec = A*D_null; % where A is the chol.
[IRF_vardec, ~, ~, ~, ~] = genIRFs(impact_vardec,0,B,0,H,sig1,sig2);
m = linspace(1,H,H);

for im = 1:length(m)     
      vardec(:,im,:) = gen_vardecomp(IRF_vardec,m(im),H);     
end







asd

% Create and Printing figures for Variance decomposition
base_path         = pwd;
which_ID          = 'vardec_GPFA_Compustat_3lags_gamgamZero_GDPDef';
print_figs        = 'no';
use_current_time  = 1; % don't save the time
plot_vardec(vardec,H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)


% Create Figure for Structural Shocks
figure
hold on
plot(Time(nlags+1:end),ssF,'LineWidth',1.5)
plot(Time(nlags+1:end),ssU,'LineWidth',1.5)
LGD = legend('Financial Shocks','Uncertainty Shocks');
LGD.FontSize = 24;
legend boxoff
axis tight
grid on