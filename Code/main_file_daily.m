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
filename                    = 'Daily';
sheet                       = 'Daily';
range                       = 'B1:L7942';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data(filename, sheet, range, do_truncation);
tf                          = isreal(dataset);
if tf == 0
      warning('Dataset has complex variables in it.')
end
dataset                     = dataset(1:end,:);
time_start                  = dataset(1,1);
time_end                    = dataset(end,1);
[~, DatasetHP]              = hpfilter(dataset(:,2:end),43200);

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

% Assess names to each variable as an array
for i = 2:size(dataset,2)
      eval([var_names{i} 'HP = DatasetHP(:,i-1);']);
end



nlags  = 24;
Y      = VXO(1+nlags:end,:);
Y2     = MoodyAaa(1+nlags:end,:);
X      = [MoodyAaa SPRV SPClose];
k      = 1;
for j = 1:nlags
      XX(:,k:k+size(X,2)-1) = [X(1+nlags-j:end-j,:)];
      k                     = k + size(X,2);
end
%XX            = [XX X(1+nlags:end,:)];% 
LM             = fitlm(XX,Y,'linear')
LM             = fitlm(XX,Y2,'linear')

[~, ~, res1]   = quick_ols(Y,XX);
[~, ~, res2]   = quick_ols(Y2,XX);

plot(Y)
hold on
plot(Y2)

corr(res1,res2)
asd






