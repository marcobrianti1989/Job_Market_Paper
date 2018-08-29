%*************************************************************************%
% Main
%
% NOTE describe variables (especially SHOCKS) in dataset
%
% last change 8/17/2018
%
% Code by Brianti, Marco e Cormun, Vito
%*************************************************************************%

clear
close all

load workspace_chol_Ufirst_GZSpread
vardec_Ulast  = vardec;
ss            = ss(:,1);

numberCPI        = strmatch('CPIInflation', system_names);
numberGDP        = strmatch('RealGDP', system_names);
numberC          = strmatch('RealCons', system_names);
numberHours      = strmatch('Hours', system_names);
numberInv        = strmatch('RealInvestment', system_names);
numberProf       = strmatch('RealProfitsaT', system_names);
numberInvent     = strmatch('RealInventories', system_names);

%numberInflation  = strmatch('Inflation', varlist);
lags             = 2;
H                = 12; %irfs horizon
mpc              = 4; %max number of principal components

% Standardize Ztilde to get one std dev shock
Ztilde  = ss;
Ztilde  = Ztilde/std(Ztilde);

% Set up the typology of transformation
logdifferences = 0;
if logdifferences == 1
      system = [nan(1,size(system,2)); diff(system)];
end
system = system(1+nlags:end,:);
pc     = pc(1+nlags:end,:);

for kk = 1:size(system,2)
      % Define inputs for local_projection
      depvarkk                    = system(:,kk);
      [~, loc_start, loc_end]     = truncate_data([depvarkk Ztilde pc]);
      loc_start                   = loc_start + lags;
      depvarkk                    = depvarkk(loc_start:loc_end);
      Ztildekk                    = Ztilde(loc_start:loc_end);
      pckk                        = pc(loc_start:loc_end,1:mpc);
      [IR{kk},~,Rsquared{kk},BL{kk},tuple{kk},VarY{kk}] = ...
            local_projection(depvarkk,pckk,Ztildekk,lags,H);
      if logdifferences == 0
            IRF(kk,:) = IR{kk};
      else
            IRF(kk,:) = cumsum(IR{kk});
      end
      % Build a table for the Variance Explained by Ztilde - Following  Stock,
      % Watson (2018) - The Economic Journal, page 928 Eq. (15)
      VarY_ih = VarY{kk};
      for ih = 1:H
            VarYY    = VarY_ih(ih);
            VarExplained(kk,ih) = sum(IRF(kk,1:ih).^2)/VarYY;
      end
      % Initiate bootstrap
      nsimul         = 2000;
      tuplekk        = tuple{kk};
      for hh = 1:H
            tuplekkhh = tuplekk{hh}; % Fix a specific horizon
            Y                             = tuplekkhh(:,1);
            X                             = tuplekkhh(:,2:end);
            [Yboot, Xboot]                = bb_bootstrap_LP(Y,X,nsimul,lags);
            for isimul = 1:nsimul
                  B                       = Xboot(:,:,isimul)'*Xboot(:,:,isimul)\...
                        (Xboot(:,:,isimul)'*Yboot(:,isimul));
                  IRF_boot(kk,hh,isimul)  = B(1);
            end
      end
end

% Select upper and lower bands
for kk = 1:size(system,2)
      IRF_bootkk = IRF_boot(kk,:,:);
      if logdifferences == 0
            IRF_boot(kk,:,:)  = IRF_bootkk;
      else
            IRF_boot(kk,:,:)  = cumsum(IRF_bootkk,2);
      end
end
IRF_boot         = sort(IRF_boot,3);
sig              = 0.05;
sig2             = 0.16;
up_bound         = floor(nsimul*sig); % the upper percentile of bootstrapped responses for CI
up_bound2        = floor(nsimul*sig2); % the upper percentile of bootstrapped responses for CI
low_bound        = ceil(nsimul*(1-sig)); % the lower percentile of bootstrapped responses for CI
low_bound2       = ceil(nsimul*(1-sig2)); % the lower percentile of bootstrapped responses for CI
IRF_up           = IRF_boot(:,:,up_bound);
IRF_up2          = IRF_boot(:,:,up_bound2);
IRF_low          = IRF_boot(:,:,low_bound);
IRF_low2         = IRF_boot(:,:,low_bound2);


%Show the graph of IRF - Figure(2)
plot2    = 1; % if plot2 = 1, figure will be displayed
n_row    = 3; % how many row in the figure
unique   = 1; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
plot_IRF_lp_unconditional(system_names,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF,H,plot2,n_row,unique)

% % Print figure authomatically if "export_figure1 = 1"
% if plot2 == 1
%       export_fig2 = 0; % if export_fig1 = 1, figure will be saved
%       export_fig_IRF_lp_unconditional(export_fig2)
% end

%Show the variance Explained - Figure(3)
plot3    = 1; % if plot2 = 1, figure will be displayed
n_row    = 3; % how many row in the figure
unique   = 1; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
plot_IRF_lp_unconditional(system_names,VarExplained,VarExplained,VarExplained,...
      VarExplained,VarExplained,H,plot3,n_row,unique)

% % Print figure authomatically if "export_figure1 = 1"
% if plot3 == 1
%       export_fig3 = 0; % if export_fig1 = 1, figure will be saved
%       export_fig_IRF_lp_unconditional(export_fig3)
% end
















