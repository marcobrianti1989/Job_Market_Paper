% -------------------------------------------------------------------------
% Generate forecast errors with autoregressive conditional mean only
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
[T,N]  = size(yt);
py     = 4;
q      = fix(4*(T/100)^(2/9));
ybetas = zeros(1+py,N);
for i = 1:N
    X    = [ones(T,1),mlags(yt(:,i),py)];
    reg  = nwest(yt(py+1:end,i),X(py+1:end,:),q);
    vyt(:,i)       = reg.resid; % forecast errors
    ybetas(:,i) = reg.beta;
end

% Save data
[T,N]  = size(vyt);
ybetas = ybetas';
dates  = 1900+(59:1/12:112-1/12)';
dates  = dates(end-T+1:end);
save arferrors dates vyt names vartype ybetas py xt

% Also write to .txt file for R code
dlmwrite('arvyt.txt',vyt,'delimiter','\t','precision',17);