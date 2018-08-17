function [B,e,V,Bs,Vs] = var_bab(Y,p,S1,S2,S3,sc,pr,bc,tr)
%% Compute bootsrapped estimates of VAR(p) for impulse responses
%  Bootsrap procedure is based on Killian (1998)
%   Input:  Y   = T-by-N matrix of data
%           p   = lag length for VAR
%           S1  = simulations for bootstrap
%           S2  = simulations for bootstrap after the bootstrap
%           S3  = simulations for bootstrap within bootstrap (optional)
%           sc  = stationarity correction, 1 = yes, 0 = no (optional)
%           pr  = print progress information, 1 = yes, 0 = no (optional)
%           bc  = bias correction before bs, 1 = yes, 0 = no (optional)
%           tr  = type of trend: const,lin,quad = 0,1,2 (optional)
%   Output: B   = OLS coefficient estimates
%           e   = OLS residual matrix
%           V   = OLS covariance matrix
%           Bs  = bootstrapped coefficient estimates
%           Vs  = bootstrapped covariance estimates
%   Author: Kyle Jurado (kej2108@columbia.edu)

%% Initialization
if nargin==4; S3=0; sc=0; pr=0; bc=1; tr=0; end;    %default: no inner loop
if nargin==5; sc=0; pr=0; bc=1; tr=0; end;          %default: no correction
if nargin==6; pr=0; bc=1; tr=0; end;                %default: no printing
if nargin==7; bc=1; tr=0; end;                      %default: bias correct
if nargin==8; tr=0; end;                            %default: constant only
gen = RandStream.create('mrg32k3a','seed',87);      %select generator
RandStream.setGlobalStream(gen);                   %set as default

%% OLS point estimation
[T,k]    = size(Y);                                 %data dimensions
if tr == 0; dt = ones(T,1); end;                    %constant only
if tr == 1; dt = [ones(T,1),(1:T)']; end;           %linear trend
if tr == 2; dt = [ones(T,1),(1:T)',(1:T)'.^2]; end; %quadratic trend
X        = [dt,mlag(Y,p)];                          %regressor matrix
X(1:p,:) = []; Y(1:p,:) = [];                       %truncate data
B        = X\Y;                                     %OLS estimate
e        = Y-X*B;                                   %OLS residuals
V        = e'*e./(T-k*p-1);                         %OLS cov matrix
I        = [eye(k*(p-1)),zeros(k*(p-1),k)];         %block diagonal
C        = [B(tr+2:end,:)';I];                      %companion matrix

%% First bootstrap
Bbar = 0;
for s = 1:S1
    ub      = e(randsample(T-p,T-p,true),:);        %draw shocks
    [Yb,Xb] = varp_gen(Y,X,B,ub,tr);                %draw bs sample
    Bbar    = Bbar + Xb\Yb./S1;                     %average estimate
end;
Phi = Bbar - B;                                     %bias estimate

%% Stationarity correction (optional)
if sc == 0; Btil = B - Phi; end;                    %no correction
if sc == 1;                                         %option for corection
m  = max(abs(eig(C)));                              %modulus of max eigval
if m >= 1; Btil = B;       end;                     %nonstationary case
if m <  1; Btil = B - Phi; end;                     %stationary case
Ctil = [Btil(tr+2:end,:)';I];                       %new companion matrix
mtil = max(abs(eig(Ctil)));                         %new modulus
if m < 1 && mtil >= 1;                              %if switched to ns case
    if pr == 1;
    disp('Stationarity correction: outer loop');    %flag for correction
    fprintf('outer iteration %d \n',s);             %display iteration
    end;
    Phic   = Phi;                                   %initialize bias
    delta  = 1;                                     %scaling term
    mtilc  = mtil;                                  %initialize new modulus
    while mtilc >= 1;                               %while nonstationary
        Phic  = delta.*Phic;                        %scale down bias
        delta = delta - 0.01;                       %smaller scale step
        Btilc = B - Phic;                           %corrected estimate
        Ctilc = [Btilc(tr+2:end,:)';I];             %companion matrix
        mtilc = max(abs(eig(Ctilc)));               %new modulus
        Btil  = Btilc;                              %new corrected estimate
    end;
end;
end;
if bc == 0; Btil = B; end;                          %no bias correction

%% Bootstrap after the Bootstrap
etil = Y - X*Btil;                                  %corrected residuals
Bs   = cell(S2,1);                                  %initialize cell array
Vs   = cell(S2,1);                                  %initialize cell array
tic;
for s = 1:S2;
    ub       = etil(randsample(T-p,T-p,true),:);    %draw shocks
    [Yb,Xb]  = varp_gen(Y,X,Btil,ub,tr);            %draw bs sample
    Bb       = Xb\Yb;                               %bs estimate
    eb       = Yb - Xb*Bb;                          %bs residuals
    Vb       = eb'*eb./(T-k*p-1);                   %bs cov matrix
    if S3 == 0; Btilb = Bb - Phi; end;              %shortcut method
    if S3 ~= 0;                                     %bs within bs (long)
    Bbar = 0;                                       %initialize mean
    for s2 = 1:S3
        ubb       = eb(randsample(T-p,T-p,true),:); %draw shocks
        [Ybb,Xbb] = varp_gen(Yb,Xb,Bb,ubb,tr);      %draw bs sample
        Bbar      = Bbar+(Xbb\Ybb)./S3;             %average estimate
    end;
    Phib = Bbar - Bb;                               %bs bias estimate
    if sc == 0; Btilb = Bb - Phib; end;             %no correction
    if sc == 1;                                     %optional correction
    m    = max(abs(eig([Bb(tr+2:end,:)';I])));      %modulus of companion
    if m >= 1; Btilb = Bb; end;                     %nonstationary case
    if m <  1; Btilb = Bb - Phib; end;              %stationary case
    mtil = max(abs(eig([Btilb(tr+2:end,:)';I])));   %new modulus
    if m < 1 && mtil >= 1;                          %if switched to ns case
        if pr == 1;
        disp('Stationarity correction: inner loop');%flag for correction
        fprintf('inner iteration %d \n',s);         %display iteration
        end;
        Phic   = Phib;                              %initialize bias
        delta  = 1;                                 %scaling term
        mtilc  = mtil;                              %initialize new modulus
        while mtilc >= 1;                           %while nonstationary
            Phic  = delta.*Phic;                    %scale down bias
            delta = delta - 0.01;                   %smaller scale step
            Btilc = Bb - Phic;                      %corrected estimate
            Ctilc = [Btilc(tr+2:end,:)';I];         %companion matrix
            mtilc = max(abs(eig(Ctilc)));           %new modulus
            Btilb = Btilc;                          %new corrected estimate
         end;
    end;
    end;
    end;
Bs{s} = Btilb;                                      %save coef sample
Vs{s} = Vb;                                         %save cov sample
if pr == 1;
if s  == 1; disp('VAR bootstrap simulations'); end;
if floor(s/round(S2/10))~=floor((s-1)/round(S2/10));
fprintf('Simulation %d, Elapsed Time %f \n',s,toc);
end
end
end
end

%% Auxiliary functions
function z = mlag(x,k,v)
% Create a matrix of n lags of a vector or matrix
%   Input:  x = matrix or vector, (nobs x k)
%           k = number of lags (default = 1)
%           v = (optional) initial values (default = 0)
%   Output: z = matrix (or vector) of lags (nobs x nvar*n)
if nargin ==1
k = 1;
v = 0;
elseif nargin == 2
v = 0;
end;
if nargin > 3
error('mlag: Wrong # of input arguments');
end;
[nobs, nvar] = size(x);
z = ones(nobs,nvar*k)*v;
for j = 1:k
    z(j+1:nobs,nvar*(j-1)+1:j*nvar) = x(1:nobs-j,:);
end
end

function [Yb,Xb] = varp_gen(Y,X,B,u,tr)
% Generate bootstrap data from resampled residuals
%   Input:   Y      = (T-p)-by-k matrix of data
%            X      = a column of ones and p lags of Y
%            B      = original OLS estimate
%            u      = resampled (T-p)-by-k residual vector
%            tr     = type of trend: const,lin,quad = 0,1,2 (optional)
%   Output: [Yb,Xb] = new generated data
%% Generate data
[Tp,k]  = size(Y);                                  %sample size
p       = (size(X,2)-(tr+1))/k;                     %number of lags in X
Y(1,:)  = X(1,:)*B+u(1,:);                          %generate first obs
for t = 2:Tp                                        %loop to generate rest
    X(t,tr+2:k+tr+1) = Y(t-1,:);                    %update first lag
    if p>1;
        X(t,tr+2+k:tr+1+p*k) = X(t-1,tr+2:tr+1+(p-1)*k);
    end                                             %update remaining lags
    Y(t,:) = X(t,:)*B+u(t,:);                       %generate next obs
end
Yb = Y; Xb = X;                                     %final output
end
