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
range                       = 'B1:T274';
do_truncation               = 1; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
dataset                     = real(dataset(1:end-1,:));

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

% Proper Transformations
percapita = 1;
if percapita == 1
      Hours       = (Hours.*Employment./Population); %Average weekly hours over population
      Consumption = (Consumption./Population);
      Investment  = (Investment./Population);
      GDP         = (GDP./Population);
      SP5001      = (SP5001./Population./GDPDef);
      SP5002      = (SP5002./Population./GDPDef);
end

% Obtaine Principal Components
filename                    = 'DatasetPC';
sheet                       = 'Quarterly';
range                       = 'C163:DD288';
dataPC                      = xlsread(filename, sheet, range);
dataPC                      = real(dataPC(1:end-1,:));
Zscore                      = 1; %remove mean and divide over the variance each series
pc                          = get_principal_components(dataPC,Zscore);
pc1                         = pc(:,1);
pc2                         = pc(:,2);
pc3                         = pc(:,3);
pc4                         = pc(:,4);

% Define the system
system_names  = {'TFPUtil','FFR','VXO','GDP','Consumption','Investment','Hours'};
for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
TFPposition = find(strcmp('TFP', system_names));
VXOposition = find(strcmp('VXO', system_names));

% Choose the correct number of lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags                    = 3;
[A,B,res,sigma]          = sr_var(system, nlags);

% Check if the VAR is stationary
test_stationarity(B');

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 2000;
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
sig1                       = 0.1;
sig2                       = 0.05;
H                          = 20;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,B,B_boot,H,sig1,sig2);

% Create and Printing figures
base_path         = pwd;
which_ID          = 'chol_';
print_figs        = 'no';
use_current_time  = 1; % don't save the time
which_shocks      = [1 2 3];
shocknames        = {'Technology Shock','Monetary Policy Shock','Uncertainty Shock'};
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
      system_names, which_ID, print_figs, use_current_time,base_path)

% Get Structural Shocks
structural_shocks = (inv(A)*res')';

% Create Figure for Structural Shocks
figure
hold on
plot(Time(nlags+1:end),structural_shocks(:,1),'LineWidth',1.5)
plot(Time(nlags+1:end),structural_shocks(:,2),'LineWidth',1.5)
plot(Time(nlags+1:end),structural_shocks(:,3),'LineWidth',1.5)
legend('Technology Shock','Monetary Policy Shock','Uncertainty Shocks')
legend boxoff
grid on

% Get variance Decomposition
[IRF_vardec, ~, ~, ~, ~] = genIRFs(A,0,B,0,H,sig1,sig2);
m = [1 4 8 16 20];
for im = 1:length(m)
      vardec(:,:,im) = gen_vardecomp(IRF_vardec,m(im),H);
end
vardec = vardec(:,1:length(which_shocks),:);




