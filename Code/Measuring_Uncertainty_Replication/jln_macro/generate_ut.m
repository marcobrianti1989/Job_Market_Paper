% -------------------------------------------------------------------------
% Compute matrix of uncertainty estimates for horizons 1 through 12
% -------------------------------------------------------------------------

% Load data
clear; clc;
load ferrors;
svf = load('svfmeans.txt');
svy = load('svymeans.txt');

% Compute objects from predictors
h   = 12;
fb  = sparse(fbetas);
thf = [svf(1,:).*(1-svf(2,:));svf(2,:);svf(3,:).^2];
xf  = svf(4:end-3,:);
gf  = svf(end-3+1:end,:);
[evf,phif] = compute_uf(xf,thf,fb,h);

% Compute uncertainty
[T,N] = size(vyt);
ut    = zeros(T,N,h);
for i = 1:N
    tic;
    yb  = sparse(ybetas(i,:));
    thy = [svy(1,i).*(1-svy(2,i));svy(2,i);svy(3,i).^2];
    xy  = svy(4:end-3,i);
    ut(:,i,:) = compute_uy(xy,thy,yb,py,evf,phif);
    fprintf('Series %d, Elapsed Time = %0.4f \n',i,toc);
end
gy = svy(end-3+1:end,:);
save ut dates ut
save geweke dates gy gf names


