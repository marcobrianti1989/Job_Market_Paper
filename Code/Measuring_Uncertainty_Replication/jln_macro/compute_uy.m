function [U] = compute_uy(xy,thy,yb,py,evf,phif)
% -------------------------------------------------------------------------
% Compute expected volatility of predictors up to horizon h
% -------------------------------------------------------------------------

% Initialize parameters
h  = length(evf);
r  = size(evf{1},2);
pf = size(phif,1)/r;
pz = (length(yb)-1-py)/r;
T  = length(xy);

% Preallocate variables
U = zeros(T,h);
if pf >1; evf0 = sparse(1,r*(pf-1),0); end;
if pf==1; evf0 = []; end;
if py >1; evy0 = sparse(1,py-1,0); end;
if py==1; evy0 = []; end;

% Construct the main phi matrix
if pf >pz; lambda_topright = sparse(1,(pf-pz)*r,0); end;
if pf==pz; lambda_topright = []; end;
if py >1;  lambda_bottom   = sparse(py-1,r*pf,0); end;
if py==1;  lambda_bottom   = []; end;
lambda   = [yb(py+2:end),lambda_topright;lambda_bottom];

phiy_top = yb(2:py+1);
if py >1; phiy_bottom = [sparse(1:py-1,1:py-1,1),sparse(py-1,1,0)]; end;
if py==1; phiy_bottom = []; end;
phiy         = [phiy_top;phiy_bottom];
phi_topright = sparse(r*pf,py,0);
phi          = [phif,phi_topright;lambda,phiy];

% Compute uncertainty using the recursion
alpha = thy(1);
beta  = thy(2);
tau2  = thy(3);
x     = xy;
for j = 1:h
    evy{j} = expectedvar(alpha,beta,tau2,x,j); %Et[(v^y_t)^2]
end;
for t = 1:T
    for j = 1:h
        ev = sparse(1:r*pf+py,1:r*pf+py,[evf{j}(t,:),evf0,evy{j}(t),evy0]);
        if j == 1; u = ev; end;
        if j  > 1; u = phi*u*phi' + ev; end; 
        U(t,j) = u(r*pf+1,r*pf+1); % select relevant entry    
    end 
end
end

% Auxiliary function
function [out] = expectedvar(a,b,t2,x,h)
% -------------------------------------------------------------------------
% Compute Et[exp{x(t+h)}] using the AR(1) law of motion for x(t)
% -------------------------------------------------------------------------
out = exp(a*(1-b^h)/(1-b)+t2/2*(1-b^(2*h))/(1-b^2)+ b^h*x);
end