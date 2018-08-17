% -------------------------------------------------------------------------
% Generate forecast errors 
% -------------------------------------------------------------------------

% Load data
clear; clc;
load jlndata; 
ind         = 132+(6:15); % "duplicate" series to remove    
data(:,ind) = []; 
names(ind)  = [];
xt          = data;

% Estimate factors
[e,fhat,lf,vf] = factors(xt,20,2,2);
[e,ghat,lg,vg] = factors(xt.^2,20,2,2);
outf     = 'mR2_fhat.out';
outg     = 'mR2_ghat.out';
[R2,mR2] = mrsq(fhat,lf,vf,names,vartype,outf);
[R2,mR2] = mrsq(ghat,lg,vg,names,vartype,outg);
ft       = [fhat,fhat(:,1).^2,ghat(:,1)]; %predictor set

% Generate forecast errors for yt
yt     = zscore(xt(:,1:132)); % only the macro data
[T,N]  = size(yt);
py     = 4;
pz     = 2;
p      = max(py,pz);
q      = fix(4*(T/100)^(2/9));
ybetas = zeros(1+py+pz*size(ft,2),N);
for i = 1:N
    X    = [ones(T,1),mlags(yt(:,i),py),mlags(ft,pz)];
    reg  = nwest(yt(p+1:end,i),X(p+1:end,:),q);
    pass = abs(reg.tstat(py+2:end)) > 2.575; % hard threshold
    keep = [ones(1,py+1)==1,pass'];
    Xnew = X(:,keep);
    reg  = nwest(yt(p+1:end,i),Xnew(p+1:end,:),q);
    vyt(:,i)       = reg.resid; % forecast errors
    ybetas(keep,i) = reg.beta;   
    fmodels(:,i)   = pass; %chosen predictors
end

% Generate AR(4) errors for ft
[T,R]  = size(ft);
pf     = 4;
q      = fix(4*(T/100)^(2/9));
fbetas = zeros(R,pf+1);
for i = 1:R
   X   = [ones(T,1),mlags(ft(:,i),pf)];
   reg = nwest(ft(pf+1:end,i),X(pf+1:end,:),q);
   vft(:,i)    = reg.resid;
   fbetas(i,:) = reg.beta';
end

% Save data
[T,N]  = size(vyt);
ybetas = ybetas';
dates  = 1900+(59:1/12:112-1/12)';
dates  = dates(end-T+1:end);
save ferrors dates vyt vft names vartype ybetas fbetas py pz pf ft xt fmodels

% Also write to .txt file for R code
dlmwrite('vyt.txt',vyt,'delimiter','\t','precision',17);
dlmwrite('vft.txt',vft,'delimiter','\t','precision',17);