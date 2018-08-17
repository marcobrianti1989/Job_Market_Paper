function z = mlags(x,k,p,v)
% -------------------------------------------------------------------------
% Create a matrix of n lags of a vector or matrix
%   Input:  x = matrix or vector, (nobs x k)
%           k = number of lags (default = 1)
%           p = number of last lags to keep (default = k)
%           v = (optional) initial values (default = 0)
%   Output: z = matrix (or vector) of lags (nobs x nvar*n)
% -------------------------------------------------------------------------
if nargin ==1 
k = 1;
v = 0;
p = k;
elseif nargin == 2
v = 0;
p = k;
elseif nargin == 3
v = 0;
end;
if p>k
error('mlags: Not enough lags'); 
end;
[nobs, nvar] = size(x);
z = ones(nobs,nvar*k)*v;
for j = 1:k
    z(j+1:nobs,nvar*(j-1)+1:j*nvar) = x(1:nobs-j,:);
end
z = z(:,end-nvar*p+1:end);
end