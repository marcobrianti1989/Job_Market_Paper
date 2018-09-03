clear
close all

filename                    = 'Monthly';
sheet                       = 'Monthly';
range                       = 'B1:AG822';
do_truncation               = 1; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
tf                          = isreal(dataset);
if tf == 0
      warning('Dataset has complex variables in it.')
end
dataset                     = real(dataset);
time_start                  = dataset(1,1);
time_end                    = dataset(end,1);

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

figure
set(gcf,'color','w');
plot(Time,zscore(VXO),'linewidth',1.5)
hold on
plot(Time,zscore(JLNUncertH1),'linewidth',2)
plot(Time,zscore(RealizedVolatility),'linewidth',1)
LEG = legend('Implied Volatility (VXO)','JLN Macro Uncertainty','Realized Volatility BDG','location','northwest');
LEG.FontSize = 24;
legend boxoff
xlabel('Years','FontSize',20)
ylabel('Standard Deviation','FontSize',20)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 16)
grid on
hold off

%invoke_export_fig('Uncertainty',[],[], cd)

cr = corr([diff(IndustrialProduction) VXO(2:end) JLNUncertH1(2:end) RealizedVolatility(2:end)]);
cr = tril(cr)
asd






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'SP500_Daily';
sheet    = 'SP500';
range    = 'A1:I17273';
[data, variables] = xlsread(filename,sheet,range);

% Assess names to each variable as an array
for i = 1:size(data,2)
      eval([variables{i} ' = data(:,i);']);
end

j = 1;
i = 2;
m = 1;
while i <= size(data,1)
      SP500_sum = 0;
      while Day(i) - Day(i-1) >= 0 && i < size(data,1)
            SP500_sum = SP500_sum + Close(i);
            i = i + 1;
      end
      SP500(j) = SP500_sum;
      monthly(j)            = m; 
      j = j + 1; 
      i = i + 1;
      if m < 12
            m = m + 1;
      else
            m = 1;
      end
end
SP500 = SP500';

xlswrite('SP_monthly',[monthly' SP500])

asdfg

clear
filename = 'SP500_Daily';
sheet    = 'SP500';
range    = 'B4:j17273';
[data, variables] = xlsread(filename,sheet,range);
day = [data(:,1)];
RVday = [data(:,end)];
j = 1;
RV = 0;
for  i = 2:length(data)
      if day(i) - day(i-1) >= 0 || i >= length(data)
            RV = RV + RVday(i);
      else           
            RVstore(j,1) = RV;
            RV = 0;
            j = j + 1;
      end
end
Monthly = monthly';
Datamonthly = [Monthly RealizedVolatility'];

xlswrite('RV_monthly2',[RVstore])

asd

j = 1;
for i = 1:length(Datamonthly)
      if floor(((i - 1)/3)) == (i - 1)/3 && length(Datamonthly) - i > 3
            RealizedVolatility_quarterly(j) = mean(RealizedVolatility(i:i+3));
            j = j + 1;
      end
end
RVq = RealizedVolatility_quarterly';

filename = 'RVq.xlsx';
xlswrite(filename,RVq)



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
load aggu

uh1 = squeeze(ut(:,:,1));
n = size(corr(uh1),1);
uncertainty = squeeze(mean(ut,2));
average_correlation = (sum(sum(abs(tril(corr(uh1))))) - n)/(n*(n-1)/2)
corr = tril(abs(corr([uncertainty uh1])));
veccorr = corr(:,1);
sorted_veccott = sort(veccorr);
Zscore = 1;
pc = get_principal_components(uh1,Zscore);

datazs = zscore(uh1);
for i = 1:size(datazs,2)
      hold on
      plot(datazs(:,i))
end

%xlswrite('uncertainty_monthly',[dates uncertainty])
%plot(dates,uncertainty(:,1))

dates(1:13)

j = 1;
for i = 1:length(uncertainty)
      if floor(((i - 1)/3)) == (i - 1)/3 && length(uncertainty) - i > 3
            u_quarterly(j,:) = mean(uncertainty(i+2:i+2+2,:),1);
            j = j + 1;
      end
end

u_quarterly = [mean(uncertainty(1:2,:),1); u_quarterly];
time = [1960.25:0.25:2011.5]';
hold on
plot(u_quarterly(:,1))

xlswrite('uncertainty',[time u_quarterly])

Inflation = [23.943
      23.917
      23.717
      23.660
      23.587
      23.767
      24.203
      24.693
      25.697
      25.947
      25.933
      26.317
      26.417
      26.487
      26.667
      26.697
      26.620
      26.720
      26.843
      26.890
      26.953
      26.910
      26.840
      26.757
      26.793
      26.757
      26.777
      26.857
      26.860
      27.037
      27.317
      27.550
      27.777
      28.013
      28.263
      28.400
      28.737
      28.930
      28.913
      28.943
      28.993
      29.043
      29.193
      29.370
      29.397
      29.573
      29.590
      29.780
      29.840
      29.830
      29.947
      29.990
      30.107
      30.220
      30.307
      30.380
      30.477
      30.533
      30.720
      30.803
      30.930
      30.980
      31.050
      31.193
      31.290
      31.490
      31.583
      31.750
      32.047
      32.337
      32.617
      32.883
      32.967
      33.167
      33.500
      33.867
      34.200
      34.533
      35.000
      35.433
      35.867
      36.433
      36.933
      37.500
      38.100
      38.633
      39.033
      39.600
      39.933
      40.300
      40.700
      41.000
      41.333
      41.600
      41.933
      42.367
      43.033
      43.933
      44.800
      45.933
      47.300
      48.567
      49.933
      51.467
      52.567
      53.200
      54.267
      55.267
      55.900
      56.400
      57.300
      58.133
      59.200
      60.233
      61.067
      61.967
      63.033
      64.467
      65.967
      67.500
      69.200
      71.400
      73.700
      76.033
      79.033
      81.700
      83.233
      85.567
      87.933
      89.767
      92.267
      93.767
      94.600
      95.967
      97.633
      97.933
      98.000
      99.133
      100.100
      101.100
      102.533
      103.500
      104.400
      105.300
      106.267
      107.233
      107.900
      109.000
      109.567
      109.033
      109.700
      110.467
      111.800
      113.067
      114.267
      115.333
      116.233
      117.567
      119.000
      120.300
      121.667
      123.633
      124.600
      125.867
      128.033
      129.300
      131.533
      133.767
      134.767
      135.567
      136.600
      137.733
      138.667
      139.733
      140.800
      142.033
      143.067
      144.100
      144.767
      145.967
      146.700
      147.533
      148.900
      149.767
      150.867
      152.100
      152.867
      153.700
      155.067
      156.400
      157.300
      158.667
      159.633
      160.000
      160.800
      161.667
      162.000
      162.533
      163.367
      164.133
      164.733
      165.967
      167.200
      168.433
      170.100
      171.433
      173.000
      174.233
      175.900
      177.133
      177.633
      177.500
      178.067
      179.467
      180.433
      181.500
      183.367
      183.067
      184.433
      185.133
      186.700
      188.167
      189.367
      191.400
      192.367
      193.667
      196.600
      198.433
      199.467
      201.267
      203.167
      202.333
      204.317
      206.631
      207.939
      210.490
      212.770
      215.538
      218.861
      213.849
      212.378
      213.507
      215.344
      217.030
      217.374
      217.297
      217.934
      219.699
      222.044
      224.568
      226.033
      227.047
      228.326
      228.808
      229.841
      231.369
      232.299
      232.045
      233.300
      234.163
      235.608
      236.839
      237.459
      236.920
      235.355
      236.912
      237.816
      237.888
      237.848
      239.452
      240.548
      242.177
      243.949
      244.010
      245.297
      247.301
      249.442
      250.468];

for ii = 1:length(Inflation)-4
      yearInflation(ii) = log(Inflation(ii+4)) - log(Inflation(ii));
end
yearInflation = [NaN; NaN; NaN; NaN; yearInflation'];

%plot(yearInflation)

filename = 'Monthly';
sheet    = 'Monthly';
range    = 'C282:C805';
[data, variables] = xlsread(filename,sheet,range);


% % Assess names to each variable as an array
% for i = 1:size(data,2)
%       eval([variables{i} ' = data(:,i);']);
% end

j = 1;
for i = 1:length(data)
      if floor(((i - 1)/3)) == (i - 1)/3 && length(data) - i > 3
            EBP_Q(j,1) = mean(data(i:i+3));
%             RealNonDurableCons_quarterly(j,1) = mean(data(i:i+3,1));
%             RealServicesCons_quarterly(j,1) = mean(data(i:i+3,2));
            j = j + 1;
      end
end

filename = 'EBP_quarterly.xlsx';

% filename = 'Cons_Decomposed.xlsx';
% xlswrite(filename,[RealServicesCons_quarterly RealNonDurableCons_quarterly RealDurableCons_quarterly])
xlswrite(filename,[EBP_Q])


















