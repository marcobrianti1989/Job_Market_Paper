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
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:AU275';
do_truncation               = 1; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
dataset                     = real(dataset(1:end-1,:));

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

% Define the system
system = [VXO RealGDP RealCons RealInvestment Hours GDPDefl M2];
system_names = {'VXO','Real GDP','Real Cons','Real Inv','Hours','GDP Defl','M2'};

% Choose the correct number of lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags           = AIC;
[A,B,res,sigma] = sr_var(system, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

% Create dataset from bootstrap
nburn             = 200;
nsimul            = 2000;
which_correction  = 'none';
blocksize         = 5;
[beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian(B,nburn,res,nsimul,which_correction,blocksize);

% Get "bootstrapped A" nsimul times
for i_simul=1:nsimul
      [A_boot(:,:,i_simul),B_boot(:,:,i_simul),~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
end

sig1                       = 0.1;
sig2                       = 0.05;
H                          = 20;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,B,B_boot,H, sig1, sig2);

%Creating and Printing figures
which_ID          = 'chol';
comment           = [which_ID '_'];
print_figs        = 'no';
use_current_time  = 1; % don't save the time
%which_shocks      = [1 2];
which_shocks      = 1;
shocknames        = {'Uncertainty Shock'};
%shocknames        = {'Technology Shock','Uncertainty Shock'};
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
      system_names, '_', print_figs, use_current_time)



