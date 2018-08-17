function [evf,phif] = compute_uf(xf,thf,fb,h)
% -------------------------------------------------------------------------
% Compute expected volatility of predictors up to horizon h,and construct
% coefficient matrix phiF
% -------------------------------------------------------------------------

% Initialize parameters
r   = size(xf,2);
pf  = size(fb,2)-1;

% Create phif matrix
phif_top = [];
for j = 2:pf+1
    phif_top = [phif_top,sparse(1:r,1:r,fb(:,j))]; 
end;
if pf > 1;
phif_bot = [sparse(1:r*(pf-1),1:r*(pf-1),1),sparse(r*(pf-1),r,0)];
phif     = [phif_top;phif_bot];
else
    phif = phif_top;
end;

% Compute evf
evf = cell(h,1);
for j = 1:h;
for i = 1:r
   alpha       = thf(1,i);
   beta        = thf(2,i);
   tau2        = thf(3,i);
   x           = xf(:,i); 
   evf{j}(:,i) = expectedvar(alpha,beta,tau2,x,j); %Et[(v^f_t)^2]
end;
end;
end

% Auxiliary function
function [out] = expectedvar(a,b,t2,x,h)
% -------------------------------------------------------------------------
% Compute Et[exp{x(t+h)}] using the AR(1) law of motion for x(t)
% -------------------------------------------------------------------------
out = exp(a*(1-b^h)/(1-b)+t2/2*(1-b^(2*h))/(1-b^2)+ b^h*x);
end