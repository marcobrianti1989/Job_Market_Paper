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
range                       = 'B1:AM274';
do_truncation               = 1; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
dataset                     = real(dataset);
time_start                  = dataset(1,1);
time_end                    = dataset(end,1);

% Obtain Principal Components
filename                    = 'DatasetPC';
sheet                       = 'Quarterly';
range                       = 'B2:DD288';
do_truncation_PC            = 1; %Do not truncate data. You will have many NaN
dataPC                      = read_data(filename,sheet,range,do_truncation_PC);
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
Zscore                      = 1; %remove mean and divide over the variance each series
pc                          = get_principal_components(dataPC(:,2:end),Zscore);
pc1                         = pc(:,1);
pc2                         = pc(:,2);
pc3                         = pc(:,3);
pc4                         = pc(:,4);

% Define the system1
system_names  = {'VXO','GDP','Consumption','Investment','Hours','TFPUtil','FFR',...
      'GovSpending','Mich5Y','Mich1Y','SP5001','GZSpread',...
      'pc1','pc2','pc3','pc4'};

for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
TFPposition = find(strcmp('TFP', system_names));
VXOposition = find(strcmp('VXO', system_names));

% Tests for lags
max_lags = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags = 2;
[A,B,res,sigma] = sr_var(system, nlags);

% Get Structural Shocks
ss = (inv(A)*res')';

y = res(:,end);
x = res(:,1:end-1);
b = (x'*y)'*(x'*x)^(-1);
sh = y - x*b';

corr(ss(:,VXOposition),res(:,VXOposition))

% Correlation between the two identifications
set(gcf,'color','white')
%area(Time(nlags+1:end),NBERDates(nlags+1:end)/10,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(Time(nlags+1:end),ss(:,VXOposition)*A(VXOposition,VXOposition),'--','linewidth',2)
plot(Time(nlags+1:end),res(:,VXOposition),'linewidth',2)
LEG = legend('VXO innovations','Constrained VXO innovations');
LEG.FontSize = 24;
legend boxoff
axis tight
grid on
hold off
close
q = 7;
k = q;
n = length(ss(:,1));
pvalue = f_test(res(:,VXOposition),ss(:,VXOposition),q,n,k);

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
which_shocks      = [1];
shocknames        = {'Uncertainty Shock'};
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
      system_names, which_ID, print_figs, use_current_time,base_path)
asd
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




