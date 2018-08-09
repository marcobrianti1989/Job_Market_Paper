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
range                       = 'B1:N274';
do_truncation               = 1; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
dataset                     = real(dataset(1:end-1,:));

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

% Proper Transformations
Hours       = (Hours.*Employment./Population); %Average weekly hours over population
Consumption = (Consumption./Population);
Investment  = (Investment./Population);
GDP         = (GDP./Population);

% Obtaining Principal Components
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
%system = [TFP VXO GDP Consumption Investment Hours GDP_Def FFR M2];
system = [TFP VXO pc1 pc2 pc3 pc4];
system_names = {'TFP','VXO','PC1','PC2','PC3','PC4'};
%system_names = {'TFP','VXO','GDP','Consumption','Investment','Hours','GDP Deflator','FFR','M2'};

% Choose the correct number of lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags           = 4;
[A,B,res,sigma] = sr_var(system, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 1000;
which_correction  = 'blocks';
blocksize         = 3;
[beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian(B,nburn,res,nsimul,which_correction,blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
      [A_boot(:,:,i_simul),B_boot(:,:,i_simul),~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
end

sig1                       = 0.1;
sig2                       = 0.05;
H                          = 40;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,B,B_boot,H,sig1,sig2);

%Creating and Printing figures
which_ID          = 'chol';
comment           = [which_ID '_'];
print_figs        = 'no';
use_current_time  = 1; % don't save the time
which_shocks      = [1 2];
%which_shocks      = 1;
%shocknames        = {'Uncertainty Shock'};
shocknames        = {'Technology Shock','Uncertainty Shock'};
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
      system_names, '_', print_figs, use_current_time)



