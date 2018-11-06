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
filename                    = 'Compustat_Data_Sep_2018_3';
sheet                       = 'WRDS';
range                       = 'A1:Q1046908';
[dataset, var_names]        = xlsread(filename, sheet, range);
tf                          = isreal(dataset);
if tf == 0
      warning('Dataset has complex variables in it.')
end

% Assess names to each variable as an array
CompanyKey    = dataset(:,1);
Quarters      = var_names(2:end,11);
Assets        = dataset(:,13);
Equity        = dataset(:,14);
CashSTInv     = dataset(:,15);
Cash          = dataset(:,16);

time = [1961:0.25:2018.25]';
year = [1961:1:2018]';
qrt  = [1:1:4]'; 
j    = 1;
for iy = 1:length(year)
      for iq = 1:length(qrt)
            index = find(strcmp(Quarters, [num2str(year(iy)) 'Q' num2str(qrt(iq))]));
            for ii = 1:length(index)
                  subCashSTInv(ii)                         = CashSTInv(index(ii));
                  subAssets(ii)                            = Assets(index(ii));
                  subEquity(ii)                            = Equity(index(ii));
                  if isnan(subCashSTInv(ii)) == 1 || isnan(subAssets(ii)) == 1 || isnan(subEquity(ii)) == 1 
                        subCashSTInv(ii)           = [];
                        subAssets(ii)              = [];
                        subEquity(ii)              = [];
                  end
            end
            % Total Assets            
            sumAssets                            = sum(subAssets);            
            AggAssets(j)                         = sumAssets;
            % Total Equity
            sumEquity                            = sum(subEquity);            
            AggEquity(j)                         = sumEquity;     
            % Cash + Short-Term Investment
            sumCashSTInv                         = sum(subCashSTInv);
            AggCashSTInv(j)                      = sumCashSTInv;
            % Counter
            j                                    = j + 1;            
      end
end

% Remove last two quarters since there are no data on that. 
AggCashSTInv        = AggCashSTInv(2:end-2);
AggEquity           = AggEquity(2:end-2);
AggAssets           = AggAssets(2:end-2);
time                = time(2:end); %this is just to show that now time is aligned with other variables

% Remove trend using the seven term henderson filter
[AggCashSTInvDt, ~] = seven_term_henderson_filter(AggCashSTInv);
[AggEquityDt, ~]    = seven_term_henderson_filter(AggEquity);
[AggAssetsDt, ~]    = seven_term_henderson_filter(AggAssets);

Cash2Assets    = (AggCashSTInvDt./AggAssetsDt)';
Cash2Equity    = (AggCashSTInvDt./AggEquityDt)';
AggEquityDt    = AggEquityDt;
avgCash2Equity = mean(Cash2Equity(24:end));
avgCash2Assets = mean(Cash2Assets(24:end));

plot(time(1:end-3),Cash2Assets(1:end-3))
hold on
plot(time(1:end-3),Cash2Equity(1:end-3))


