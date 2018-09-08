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
range                       = 'B1:AT274';
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

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} 'HP1S = DatasetHP1S(:,i);']);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HP1S w/o EBP
Y         = CashFlowHP1S;
X         = [MacroUncertH1HP1S GDPHP1S];
LM        = fitlm(X,Y,'linear')
% HP1S with EBP
Y         = CashFlowHP1S;
X         = [MacroUncertH1HP1S EBPHP1S GDPHP1S];
LM        = fitlm(X,Y,'linear')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HP w/o EBP
Y         = CashFlowHP;
X         = [MacroUncertH1HP GDPHP];
LM        = fitlm(X,Y,'linear')
% HP with EBP
Y         = CashFlowHP;
X         = [MacroUncertH1HP EBPHP GDPHP];
LM        = fitlm(X,Y,'linear')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlags = 6;
% diff w/o EBP
clear X
Y = CashFlow(1+nlags:end);
k = 1;
for j = 1:nlags+1
      X(:,k:k+1) = [MacroUncertH1(1+nlags-j+1:end-j+1)...
            GDP(1+nlags-j+1:end-j+1)];
      k = k + 2;
end
LM        = fitlm(X,Y,'linear')
% diff with EBP
clear X
Y = CashFlow(1+nlags:end);
k = 1;
for j = 1:nlags+1
      X(:,k:k+2) = [MacroUncertH1(1+nlags-j+1:end-j+1)...
            GDP(1+nlags-j+1:end-j+1) EBP(1+nlags-j+1:end-j+1)];
      k = k + 3;
end
LM        = fitlm(X,Y,'linear')





