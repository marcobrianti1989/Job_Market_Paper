% -------------------------------------------------------------------------
% Generate forecast errors with no predictors
% -------------------------------------------------------------------------

% Load data
clear; clc;
load jlndata; 
ind         = 132+(6:15); % "duplicate" series to remove    
data(:,ind) = []; 
names(ind)  = [];
xt          = data;

% Generate forecast errors for yt
yt     = zscore(xt(:,1:132)); % only the macro data
vyt    = yt;

% Save data
[T,N]  = size(vyt);
dates  = 1900+(59:1/12:112-1/12)';
dates  = dates(end-T+1:end);
save npferrors dates vyt names vartype xt

% Also write to .txt file for R code
dlmwrite('npvyt.txt',vyt,'delimiter','\t','precision',17);