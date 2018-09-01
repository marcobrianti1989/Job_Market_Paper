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
filename                    = 'Quarterly';
sheet                       = 'Quarterly Data';
range                       = 'B1:AR274';
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

% Obtain Principal Components
filename                    = 'DatasetPC';
sheet                       = 'Quarterly';
range                       = 'B2:DA288';
do_truncation_PC            = 1; %Do not truncate data. You will have many NaN
dataPC                      = read_data(filename,sheet,range,do_truncation_PC);
tfPC                        = isreal(dataset);
if tfPC == 0
      warning('DatasetPC has complex variables in it.')
end
dataPC                      = real(dataPC);
time_start_PC               = dataPC(1,1);
time_end_PC                 = dataPC(end,1);

% Align the two datasets
align_datasets = 1;
if align_datasets == 1
      if time_start < time_start_PC
            loc_start = find(dataset(:,1) == time_start_PC);
            dataset = dataset(loc_start:end,:);
      elseif time_start > time_start_PC
            loc_start = find(dataPC(:,1) == time_start);
            dataPC = dataPC(loc_start:end,:);
      end
      if time_end < time_end_PC
            loc_end = find(dataPC(:,1) == time_end);
            dataPC = dataPC(1:loc_end,:);
      elseif time_end > time_end_PC
            loc_end = find(dataset(:,1) == time_end);
            dataset = dataset(1:loc_end,:);
      end
end

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
      Hours         = Hours + Employment - Population; %Average weekly hours over population
      Consumption   = NonDurableCons + ServiceCons - Population;
      Investment    = Investment + DurableCons - Population;
      GDP           = GDP - Population;
      SP5001        = SP5001 - Population - GDPDef;
      SP5002        = SP5002 - Population - GDPDef;
      GovPurchases  = GovPurchases - Population;
      CashFlow      = CashFlow - Population - GDPDef;
end

% Obtaine Principal Components
Zscore                      = 1; %remove mean and divide over the variance each series
pc                          = get_principal_components(dataPC(:,2:end),Zscore);
pc1                         = pc(:,1);
pc2                         = pc(:,2);
pc3                         = pc(:,3);
pc4                         = pc(:,4);

% Define the system1
% system_names  = {'SP5001','MacroUncertH1','TFPUtil','GDP','Consumption',...
%       'Investment','Hours','YearInflation','FFR','GovSpending','CapUtilization','Inventories'};
system_names  = {'MacroUncertH1','GDP','Consumption','Investment',...
      'Hours','EBP','SP5002'};

for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
TFPposition = find(strcmp('TFP', system_names));
VXOposition = find(strcmp('VXO', system_names));
Uposition   = find(strcmp('MacroUncertH1', system_names));

% Tests for lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags           = 4;
[A,B,res,sigma] = sr_var(system, nlags);

% Get Structural Shocks
ss = (inv(A)*res')';

% Test for sufficient information - H0: regressors on PC are equal to zero.
reg_PC                     = [pc1 pc2 pc3 pc4];
reg_PC                     = reg_PC(1+nlags:end,:);
ushock_restricted          = res(:,2);
k                          = size(reg_PC,2);
q                          = size(reg_PC,2);
[~,~,ushock_unrestricted]  = quick_ols(ushock_restricted,reg_PC);
TT                         = size(reg_PC,1);
pvalue_FGtest              = f_test(ushock_restricted,ushock_unrestricted,q,TT,k);

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
H                          = 20;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,B,B_boot,H,sig1,sig2);

% Create and Printing figures
base_path         = pwd;
which_ID          = 'chol_';
print_figs        = 'no';
use_current_time  = 1; % don't save the time
which_shocks      = [Uposition];
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







