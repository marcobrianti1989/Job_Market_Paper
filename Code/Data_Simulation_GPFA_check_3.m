clear
close all

% x  = rho*x(-1) + rh0*y(-1) + rho*z(-2) + epsU + epsF;
% y  = rho*x(-1) + rh0*y(-1) + rho*z(-2) + epsU + epsF;
% z  = rho*x(-1) + rh0*y(-1) + rho*z(-2) + epsU - epsF;

% Parameterization
rhox  = 0.8;
rhoy  = 0.6;
rho   = 0.2;
rhoz  = 0.4;
rhozz = 0.6;
sigU  = 1;
sigF  = 1;
sigC  = 1;
A11   = 1;
A12   = 3.5;
A21   = 0.75;
A22   = 0.5;

% Number of Observations
T     = 100000;

% Predetermined matrices
U = zeros(T,1);
F = zeros(T,1);
C = zeros(T,1);
shocksU = sigU*randn(T,1);
shocksF = sigF*randn(T,1);
shocksC = sigC*randn(T,1);

% Structural Economy
% for i = 2:T
%       U(i) = rhox*U(i-1) + rho*F(i-1) + 2*shocksU(i) + shocksF(i);
%       F(i) = rho*U(i-1) + rhoy*F(i-1) + shocksU(i) + 2*shocksF(i) + shocksC(i);
%       C(i) = rho*U(i-1) - rho*F(i-1)  + rhoz*C(i-1) ...
%             + shocksU(i) - shocksF(i) + shocksC(i);
% end

for i = 2:T
      U(i) = rhox*U(i-1) + rhoz*C(i-1)   + A11*shocksU(i)  + A21*shocksF(i);  
      C(i) = rho*U(i-1)  + rhozz*C(i-1)  + A21*shocksU(i)  - A22*shocksF(i);
end

% Define the system1
%system_names  = {'U','F','C'};
system_names  = {'U','C'};

for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
Uposition   = find(strcmp('U', system_names));
EBPposition = find(strcmp('F', system_names));
IVposition  = find(strcmp('C', system_names));

% Cholesky decomposition
nlags           = 1;
[A,B,res,sigma] = sr_var(system, nlags);

% Generalized Penalty Function Approach
objgam     = 1;
delta      = 0;
threshold  = 0.1;
add        = 0.01;
while objgam >= 10^(-16)
      while objgam >= threshold
            warning off
            SRhorizon       = 1;
            SRhorizonIV     = 4;
            [impact, gamma] = identification_GPFA(A,B,SRhorizon,...
                  SRhorizonIV,Uposition,...
                  Uposition,IVposition,-delta);
            %       [impact, gamma] = identification_GPFA_2delta(A,B,SRhorizon,SRhorizonIV,...
            %             EBPposition,Uposition,IVposition,delta1,delta);
            gamF   = gamma(:,1);
            gamU   = gamma(:,2);
            objgam = gamF'*gamU
            delta  = delta + add;
      end
      threshold = threshold/10;
      add = add/10;
end
delta = -delta; % Negative because the first variable is EBP and cash respons negatively on impact

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 10;
which_correction  = 'none';
blocksize         = 4;
[beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian(B,nburn,res,nsimul,which_correction,blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
      % Cholesky decomposition
      [A_boot(:,:,i_simul),B_boot(:,:,i_simul),~,~] = ...
            sr_var(data_boot2(:,:,i_simul), nlags);
      % GPFA identification strategy
      warning off
      [impact_boot(:,:,i_simul), gamma_boot(:,:,i_simul)] = ...
            identification_GPFA(A_boot(:,:,i_simul),B_boot(:,:,i_simul),...
            SRhorizon,SRhorizonIV,Uposition,Uposition,IVposition,delta);
      i_simul
end

% Generate IRFs with upper and lower bounds
sig1                       = 0.05;
sig2                       = 0.025;
H                          = 20;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(impact,impact_boot,B,B_boot,H,sig1,sig2);

% Create and Printing figures for IRFs
base_path         = pwd;
which_ID          = 'GPFA_Compustat_3lags_gamgamZero_GDPDef';
print_figs        = 'no';
use_current_time  = 1; % (don't) save the time
which_shocks      = [1 2]; %[Uposition];
shocknames        = {'Financial Shock','Uncertainty Shock'};
plot_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)

% Get variance Decomposition
N = null(gamma');
D_null = [gamma N];
impact_vardec = A*D_null; % where A is the chol.
[IRF_vardec, ~, ~, ~, ~] = genIRFs(impact_vardec,0,B,0,H,sig1,sig2);
m = linspace(1,H,H);
for im = 1:length(m)
      vardec(:,im,:) = gen_vardecomp(IRF_vardec,m(im),H);
end

% Theoretical IRF - U shock
xtheory = zeros(H,2);
ytheory = zeros(H,2);

xtheory(1,2) = 2*sigU;
ytheory(1,2) = sigU;
ztheory(1,2) = sigU;
for i = 2:H
      xtheory(i,2) = rhox*xtheory(i-1,2) + rho*ytheory(i-1,2);
      ytheory(i,2) = rho*xtheory(i-1,2) + rhoy*ytheory(i-1,2);
      ztheory(i,2) = rho*xtheory(i-1,2) - rho*ytheory(i-1,2)  + rhoz*ztheory(i-1,2);
end
% Theoretical IRF - F shock
xtheory(1,1) = sigF;
ytheory(1,1) = 2*sigF;
ztheory(1,1) = -sigF;
for i = 2:H
      xtheory(i,1) = rhox*xtheory(i-1,1) + rho*ytheory(i-1,1);
      ytheory(i,1) = rho*xtheory(i-1,1) + rhoy*ytheory(i-1,1);
      ztheory(i,1) = rho*xtheory(i-1,1) - rho*ytheory(i-1,1)  + rhoz*ztheory(i-1,1);
end

theoryIRF = zeros(3,H,2);
theoryIRF(1,:,:) = xtheory;
theoryIRF(2,:,:) = ytheory;
theoryIRF(3,:,:) = ztheory;

print_figs = 'no';
plot_IRFs_Empirical_Theoretical_2CIs(IRFs,ub1,lb1,ub2,lb2,theoryIRF,H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)
print_figs = 'no';
asd
% Create and Printing figures for Variance decomposition
base_path         = pwd;
which_ID          = 'vardec_GPFA_Compustat_3lags_gamgamZero_GDPDef';
print_figs        = 'no';
use_current_time  = 1; % don't save the time
plot_vardec(vardec,H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)

% Get Structural Shocks
ss = get_structural_shocks_general(A,gamma,res,which_shocks);
ssF = ss(:,1);
ssU = ss(:,2);
corr(ssF,shocksF(2:end))
corr(ssU,shocksU(2:end))

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
