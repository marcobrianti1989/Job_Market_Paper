clear
close all

% Create the correct path
base_path = pwd;
if exist([base_path '\Data'], 'dir')
      addpath([base_path '\Data']) %for Microsoft
else
      addpath([base_path '/Data']) %for Mac
end

% Reading Data
filename                    = 'Compustat_Data_Sep_2018_2';
sheet                       = 'WRDS';
range                       = 'A1:Q1046908';
[dataset, var_names]        = xlsread(filename, sheet, range);
tf                          = isreal(dataset);
if tf == 0
      warning('Dataset has complex variables in it.')
end

% Assess names to each variable as an array
CompanyKey    = dataset(:,1);
Quarters      = var_names(2:end,10);
Assets        = dataset(:,12);
Equity        = dataset(:,13);
CashSTInv     = dataset(:,14);
Cash          = dataset(:,15);



time = [1961:0.25:2018.25]';
year = [1961:1:2018]';
qrt  = [1:1:4]';
j    = 1;
for iy = 1:length(year)
      for iq = 1:length(qrt)
            index = find(strcmp(Quarters, [num2str(year(iy)) 'Q' num2str(qrt(iq))]));
            % Cash
            subCash                              = Cash(index);
            subCash(isnan(subCash))              = [];            
            sumCash                              = sum(subCash);
            AggCash(j)                           = sumCash;
            % Total Assets          
            subAssets                            = Assets(index);
            subAssets(isnan(subAssets))          = [];            
            sumAssets                            = sum(subAssets);            
            AggAssets(j)                         = sumAssets;
            % Total Equity
            subEquity                            = Equity(index);
            subEquity(isnan(subEquity))          = [];
            sumEquity                            = sum(subEquity);            
            AggEquity(j)                         = sumEquity;     
            % Cash + Short-Term Investment
            subCashSTInv                         = CashSTInv(index);
            subCashSTInv(isnan(subCashSTInv))    = [];
            sumCashSTInv                         = sum(subCashSTInv);
            AggCashSTInv(j)                      = sumCashSTInv;
            % Counter
            j                                    = j + 1;
            
      end
end
AggCash        = AggCash(1:end-2);
AggCashSTInv   = AggCashSTInv(1:end-2);
AggEquity      = AggEquity(1:end-2);
AggAssets      = AggAssets(1:end-2);

[AggCashSTInvDt, ~] = seven_term_henderson_filter(AggCashSTInv);
[AggCashDt, ~]      = seven_term_henderson_filter(AggCash);
[AggEquityDt, ~]    = seven_term_henderson_filter(AggEquity);
[AggAssetsDt, ~]    = seven_term_henderson_filter(AggAssets);

NCash = AggCashSTInvDt./AggAssetsDt;

plot(time(1:end-3),NCash(1:end-3))

