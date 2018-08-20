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
range                       = 'B1:AK274';
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

% Define the system1
system_names1  = {'VXO','GDP','Consumption','Investment','Hours','TFPUtil','FFR','GovSpending','pc1','pc2','pc3','pc4'};

for i = 1:length(system_names1)
      system1(:,i) = eval(system_names1{i});
end
TFPposition1 = find(strcmp('TFP', system_names1));
VXOposition1 = find(strcmp('VXO', system_names1));

% Choose the correct number of lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system1,max_lags);

% Cholesky decomposition
nlags                        = 3;
[A1,B1,res1,sigma1]          = sr_var(system1, nlags);

% Check if the VAR is stationary
test_stationarity(B1');

% Get Structural Shocks
ss1 = (inv(A1)*res1')';

% regreg = (res1(:,2)'*res1(:,1))/(res1(:,1)'*res1(:,1));
% ss1_con = res1(:,2)/A1(2,2) - res1(:,1)*regreg/A1(2,2);
% q = 7;
% k = q;
% n = length(ss1(:,1));
% pvalue = f_test(ss1_con,ss1(:,1), q, n, k);

% Define the system2
system_names2  = {'GDP','Consumption','Investment','Hours','TFPUtil','FFR','GovSpending','VXO','pc1','pc2','pc3','pc4'};
%system_names2  = {'GDP','VXO'};

for i = 1:length(system_names2)
      system2(:,i) = eval(system_names2{i});
end
TFPposition2 = find(strcmp('TFP', system_names2));
VXOposition2 = find(strcmp('VXO', system_names2));

% Cholesky decomposition
[A2,B2,res2,sigma2] = sr_var(system2, nlags);

% Get Structural Shocks
ss2 = (inv(A2)*res2')';

y = res2(:,end);
x = res2(:,1:end-1);

b = (x'*y)'*(x'*x)^(-1);

sh = y - x*b';

corr(ss2(:,end-4),res1(:,1))

% Correlation between the two identifications
corr(res1(:,1),ss2(:,end))
plot(Time(4:end),res1(:,1)/A1(1,1),'linewidth',2)
hold on
plot(Time(4:end),ss2(:,end-4),'linewidth',2)
legend('VXO first','VXO last')

q = 7;
k = q;
n = length(ss1(:,1));
pvalue = f_test(y,sh,q,n,k);

asdf

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




