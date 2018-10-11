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
range                       = 'B1:BG274';
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
for iif = 3:length(dataset)
      [~ , dataset_HPgen]   = hpfilter(dataset(1:iif,:),1600);
      DatasetHP1S(iif,:)    = dataset_HPgen(end,:);
end

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} 'HP = DatasetHP(:,i);']);
end

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} 'HP1S = DatasetHP1S(:,i);']);
end

% Proper Transformations - All the variables should be in logs
percapita = 1;
if percapita == 1
      Hours               = Hours + Employment - Population; %Average weekly hours over population
      Consumption         = NonDurableCons + ServiceCons - Population;
      Investment          = Investment + DurableCons - Population;
      GDP                 = GDP - Population;
      SP5001              = SP5001 - Population - GDPDef;
      SP5002              = SP5002 - Population - GDPDef;
      GovPurchases        = GovPurchases - Population;
      CPAdj               = CorpProfitsAdj - Population - GDPDef;
      TenYTreasury        = TenYTreasury;
      SecuAssets          = CorpSecuritiesAssets - Population - GDPDef;
      M2                  = M2 - Population;
else
      Consumption         = NonDurableCons + ServiceCons;
      Investment          = Investment + DurableCons;
      SP5001              = SP5001 - GDPDef;
      SP5002              = SP5002 - GDPDef;
      CPAdj               = CorpProfitsAdj - GDPDef;
      SecuAssets          = CorpSecuritiesAssets - GDPDef;
end
% Other Transformations
CashFlow            = CashFlow - CorpProfitsAdj;
Dividends           = Dividends - CorpProfitsAdj;
ProfTransfers       = ProfitsTransfers - CorpProfitsAdj;
ConsFixedK          = ConsFixedK - CorpProfitsAdj;
Cash                = Cash - GDPDef;

% Define the system1
system_names  = {'EBP','MacroUncertH3','GDP','Consumption',...
      'Investment','Hours','Cash2Assets','GDPDef'};

for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
Uposition   = find(strcmp('MacroUncertH3', system_names));
EBPposition = find(strcmp('EBP', system_names));
IVposition  = find(strcmp('Cash2Assets', system_names));
GDPposition = find(strcmp('GDP', system_names));
Cposition   = find(strcmp('Consumption', system_names));

% Tests for lags
max_lags     = 4;
[AIC,BIC,HQ] = aic_bic_hq(system,max_lags);

% Cholesky decomposition
nlags           = 3;
[A,B,res,sigma] = sr_var(system, nlags);

% Generalized Penalty Function Approach
objgam = 1;
delta = 4.7121851891011;
while abs(objgam) >= 0.05;%10^(-8)
      warning off
      SRhorizon       = 1;
      SRhorizonIV     = 4;
      [impact, gamma] = identification_GPFA(A,B,SRhorizon,SRhorizonIV,EBPposition,...
            Uposition,IVposition,-delta);
      gamF   = gamma(:,1);
      gamU   = gamma(:,2);
      objgam = gamF'*gamU
      delta  = delta + 0.000000000000001
end
delta = - delta; % Negative because the first variable is EBP and cash respons negatively on impact

% Create dataset from bootstrap
nburn             = 0;
nsimul            = 200;
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
            SRhorizon,SRhorizonIV,EBPposition,Uposition,IVposition,delta);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read Ramey Shocks
filename                    = 'RameyShocks';
sheet                       = 'Sheet1';
range                       = 'B1:M278';
do_truncation_RameyShocks   = 0; %Do not truncate data. You will have many NaN
[Rameyshocks, Rameynames]   = read_data(filename,sheet,range,do_truncation_RameyShocks);
tfRameyShocks               = isreal(Rameyshocks);
if tfRameyShocks == 0
      warning('DatasetPC has complex variables in it.')
end

% Assess names to each shock as an array
for i = 1:size(Rameyshocks,2)
      eval([Rameynames{i} ' = Rameyshocks(:,i);']);
end

for i = 2:size(Rameyshocks,2)
      % Get Structural Shocks
      ss = get_structural_shocks_general(A,gamma,res,which_shocks);
      ssF = ss(:,1);
      ssU = ss(:,2);
      [tuple, ~, ~] = truncate_data([TimeShocks Rameyshocks(:,i)]);
      % Align the two datasets
      time_start_Ramey = tuple(1,1);
      time_end_Ramey   = tuple(end,1);
      if Time(nlags+1) < time_start_Ramey && Time(end) > time_end_Ramey
            disp([num2str(i),' is Case 1'])
            loc_start_Ramey = find(Time(:) == time_start_Ramey);
            loc_end_Ramey   = find(Time(:) == time_end_Ramey);
            ssU             = ssU(loc_start_Ramey:loc_end_Ramey);
            ssF             = ssF(loc_start_Ramey:loc_end_Ramey);
      elseif Time(nlags+1) > time_start_Ramey && Time(end) < time_end_Ramey
            disp([num2str(i),' is Case 2'])
            loc_start_Ramey = find(tuple(:,1) == Time(nlags+1));
            loc_end_Ramey   = find(tuple(:,1) == Time(end));
            tuple           = tuple(loc_start_Ramey:loc_end_Ramey,:);
      elseif Time(nlags+1) < time_start_Ramey && Time(end) < time_end_Ramey
            disp([num2str(i),' is Case 3'])
            loc_start_Ramey = find(Time(:) == time_start_Ramey);
            ssU             = ssU(loc_start_Ramey:end);
            ssF             = ssF(loc_start_Ramey:end);            
            loc_end_Ramey   = find(tuple(:,1) == Time(end));
            tuple           = tuple(1:loc_end_Ramey,:);
      elseif Time(nlags+1) > time_start_Ramey && Time(end) > time_end_Ramey
            disp([num2str(i),' is Case 4'])
            loc_start_Ramey = find(tuple(:,1) == Time(nlags+1));
            tuple           = tuple(loc_start_Ramey:end,:);
            loc_end_Ramey   = find(Time(:) == time_end_Ramey);
            ssU             = ssU(1+nlags:loc_end_Ramey);
            ssF             = ssF(1+nlags:loc_end_Ramey);
      end
      bigtuple                = [tuple ssU ssF];
      Rameynames{i}
      [corri, pvaluei]        = corrcoef(bigtuple(:,2:end))
      cell_tuple              = num2cell(bigtuple);
      cell_tuple_store{i-1}   = cell_tuple;
      
end

