function [ehat,Fhat,lamhat,ve2] = factors(X,kmax,jj,DEMEAN)
% -------------------------------------------------------------------------
% Estimate latent factors from observable data X using PCA
%   Input:  X      = observable data
%           kmax   = maximum number of factors to consider
%           jj     = information criterion from Bai and Ng (2002)
%           DEMEAN = 0 - no, 1 - yes, 2 - standardize
%   Output: ehat   = idiosyncratic errors
%           Fhat   = latent factor estimates
%           lamhat = factor loadings
%           ve2    = eigenvalues of data covariance matrix
% -------------------------------------------------------------------------

[ic1,chat,fhat,eigval]  = nbplog(X,kmax,jj,DEMEAN); 
icstar    = ic1;
R2_static = sum(eigval(1:icstar))/sum(eigval);
if DEMEAN == 2;
[ehat,Fhat,lamhat,ve2]  = pc(standard(X),icstar);
end
if DEMEAN == 1;
[ehat,Fhat,lamhat,ve2]  = pc(X-repmat(mean(X),size(X,1),1),icstar); 
end
if DEMEAN == 0;
[ehat,Fhat,lamhat,ve2]  = pc(X,icstar); 
end
end

% Auxiliary functions
function [ic1, chat,Fhat,eigval]=nbplog(x,kmax,jj,DEMEAN);
%function [ic1, lambda,Fhat,IC1]=nbplog(x,kmax,jj,DEMEAN);
T=size(x,1);
N=size(x,2);
NT=N*T;
NT1=N+T;
CT=zeros(1,kmax);
ii=1:1:kmax;
if jj ==1;CT(1,:)=log(NT/NT1)*ii*NT1/NT;end;
if jj==2; CT(1,:)=(NT1/NT)*log(min([N;T]))*ii;end;
GCT=min([N;T]);
if jj==3; CT(1,:)=ii*log(GCT)/GCT; end;
if jj==4; CT(1,:)=2*ii/T; end;
if jj==5; CT(1,:)=log(T)*ii/T;end;
if jj==6; CT(1,:)=2*ii*NT1/NT; end;
if jj==7; CT(1,:)=log(NT)*ii*NT1/NT;end;
if jj==8; CT(1,:)= 2*ii*(sqrt(N)+sqrt(T))^2/(NT); end;% new modified CP
  
  
 if DEMEAN ==2;
 X=standard(x);
 end;

if DEMEAN ==1;
 X=x-repmat(mean(x),T,1);
 end;
if DEMEAN==0;
  X=x;;
  end;
IC1=zeros(size(CT,1),kmax+1);
Sigma=zeros(1,kmax+1);
if T< N;
[ev,eigval,ev1]=svd(X*X');
sumeigval=cumsum(diag(eigval))/sum(diag(eigval));
Fhat0=sqrt(T)*ev;
Lambda0=X'*Fhat0/T;
else;
[ev,eigval,ev1]=svd(X'*X);
sumeigval=cumsum(diag(eigval))/sum(diag(eigval));
Lambda0=sqrt(N)*ev;
Fhat0=X*Lambda0/N;
end;



if jj <= 8;
for i=kmax:-1:1;
Fhat=Fhat0(:,1:i);
%lambda=Fhat'*X;
lambda=Lambda0(:,1:i);
chat=Fhat*lambda';
%disp([i sumeigval(i) sum(sum(chat.*chat))/sum(sum(X.*X))]);
ehat=X-chat;
Sigma(i)=mean(sum(ehat.*ehat/T));
IC1(:,i)=log(Sigma(i))+CT(:,i);
end;
Sigma(kmax+1)=mean(sum(X.*X/T));
IC1(:,kmax+1)=log(Sigma(kmax+1));
ic1=minindc(IC1')';
ic1=ic1 .*(ic1 <= kmax);
end;
if jj==9;

  for j=1:rows(sumeigval);
    if sumeigval(j) >= .5; ic1=j; break; end;
  end;    
 end; 
Fhat=[];
Fhat=Fhat0(:,1:kmax);
Lambda=Lambda0(:,1:kmax);
chat=Fhat*Lambda';
eigval=diag(eigval);
end

function [ehat,fhat,lambda,ss]=pc(y,nfac);

[bigt,bign]=size(y);
yy=y'*y;
[Fhat0,eigval,Fhat1]=svd(yy);
lambda=Fhat0(:,1:nfac)*sqrt(bign);
fhat=y*lambda/bign;
ehat=y-fhat*lambda';

ve2=sum(ehat'.*ehat')'/bign;
ss=diag(eigval);
end

function pos=minindc(x);
ncols=size(x,2);
nrows=size(x,1);
pos=zeros(ncols,1);
seq=seqa(1,1,nrows);
for i=1:ncols;
dum=min(x(:,i));
dum1= seq .* ( (x(:,i)-dum) ==0);
pos(i)=sum(dum1);
end;
end

function x=standard(y);
T=size(y,1);
N=size(y,2);
my=repmat(mean(y),T,1);
sy=repmat(std(y),T,1);
x=(y-my)./sy;

%x=(y-kron(mean(y),ones(rows(y),1)))./kron(std(y),ones(rows(y),1));
end

function seq=seqa(a,b,c);
seq=(a:b:(a+b*(c-1)))';
end