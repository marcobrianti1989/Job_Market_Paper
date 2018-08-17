% -------------------------------------------------------------------------
% Compute matrix of uncertainty estimates (ar conditional mean only)
% -------------------------------------------------------------------------

% Load data
clear; clc;
load arferrors;
svy = load('arsvmeans.txt');

% Compute uncertainty
h = 12;
[T,N] = size(vyt);
thy   = [svy(1,:).*(1-svy(2,:));svy(2,:);svy(3,:).^2];
xy    = svy(4:end-3,:);
gy    = svy(end-3+1:end,:);
ut    = zeros(T,N,h);
py    = size(ybetas,2)-1;
for i = 1:N
    tic;
    phi   = sparse([ybetas(i,2:end);[eye(py-1),zeros(py-1,1)]]);
    a     = thy(1,i);
    b     = thy(2,i);
    t2    = thy(3,i);
    x     = xy(:,i);
    for j = 1:h
        evy{j} =  exp(a*(1-b^j)/(1-b)+t2/2*(1-b^(2*j))/(1-b^2)+ b^j*x);
    end
    for t = 1:T
        for j = 1:h
            ev = sparse(1,1,evy{j}(t),py,py);
            if j == 1; u = ev; end;
            if j  > 1; u = phi*u*phi' + ev; end;
            ut(t,i,j) = u(1,1);
        end
    end
    fprintf('Series %d, Elapsed Time = %0.4f \n',i,toc);
end

% Save results
save arut dates ut
save argeweke dates gy names