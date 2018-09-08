%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Marco Brianti, PhD Candidate, Boston College, Department of Economics, August 8, 2018
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% Reading Data
filename                    = 'Monthly';
sheet                       = 'Monthly';
range                       = 'B1:AG828';
do_truncation               = 1; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
tf                          = isreal(dataset);
if tf == 0
      warning('Dataset has complex variables in it.')
end
dataset                     = real(dataset);
time_start                  = dataset(1,1);
time_end                    = dataset(end,1);
[~, DatasetHP]              = hpfilter(dataset,1600);

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} 'HP = DatasetHP(:,i);']);
end

% Proper Transformations - All the variables should be in logs
percapita = 1;
if percapita == 1
      Consumption   = RealNonDurableCons + RealServicesCons - Population;
      Investment    = RealDurableCons - Population;
      IP            = IndustrialProduction - Population;
      SP500         = SP500 - Population - CPI;
end

% Define the system1
% system_names  = {'SP5001','MacroUncertH1','TFPUtil','GDP','Consumption',...
%       'Investment','Hours','YearInflation','FFR','GovSpending','CapUtilization','Inventories'};
system_names  = {'EBP','IP','Consumption','Investment','SP500','UnRate'};

for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
Uposition   = find(strcmp('JLNUncertH1', system_names));
EBPposition = find(strcmp('EBP', system_names));

% Tests for lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags           = 4;
[A,B,res,sigma] = sr_var(system, nlags);

j = 0;
resIP = res(1:end-j,2);
resC  = res(1:end-j,3);
resSP = res(1:end-j,6);

% Linear Regression Model. It is the same of regress but it provides much
% more useful information.
Y = JLNUncertH1(1+nlags+j:end);
Y2 = VXOHP(1+nlags+j:end);
IPGrowth = diff(IP);
%autocorr(IPGrowth.^2)
IPGrowth = IPGrowth(nlags:end-j);
X = [resIP IPGrowth];
%LM = fitlm(X,Y,'linear')

% Get Structural Shocks
ss = (inv(A)*res')';

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 500;
which_correction  = 'none';
blocksize         = 4;
[beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian(B,nburn,res,nsimul,which_correction,blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
      % Cholesky decomposition
      [A_boot(:,:,i_simul),B_boot(:,:,i_simul),~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
end

% Generate IRFs with upper and lower bounds
sig1                       = 0.05;
sig2                       = 0.025;
H                          = 120;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,B,B_boot,H,sig1,sig2);

% Create and Printing figures
base_path         = pwd;
which_ID          = 'chol_';
print_figs        = 'no';
use_current_time  = 1; % don't save the time
which_shocks      = [1];
shocknames        = {'Uncertainty Shock'};
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)

% Get variance Decomposition
[IRF_vardec, ~, ~, ~, ~] = genIRFs(A,0,B,0,H,sig1,sig2);
m = linspace(1,H,H);
for im = 1:length(m)
      vardec(:,:,im) = gen_vardecomp(IRF_vardec,m(im),H);
end
vardec = vardec(:,which_shocks,:);
horz = linspace(0,H,H);
figure
hold on
plot(horz,vardec(2,:),'linewidth',2)
grid on
legend boxoff
xlabel('Horizon')
ylabel('Variance Explained')
title('Variance Explained Of Real GDP')

%save workspace_nicespecification_cons_inv_adjusted_Ulast







