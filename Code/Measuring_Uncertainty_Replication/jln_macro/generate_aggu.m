% -------------------------------------------------------------------------
% Generate aggregate uncertainty estimates
% -------------------------------------------------------------------------

% Initialization
clear; clc;
load ut;
[T,N,h] = size(ut);

% Cross-sectional average
utcsa = squeeze(mean(sqrt(ut),2));

% Principal component analysis
utpca = zeros(T,h);
for j = 1:h
   logu    = log(sqrt(ut(:,:,j)));
   dlogu(:,:,j)  = logu(2:end,:) - logu(1:end-1,:);
   [de,du,dl,dv] = factors(dlogu(:,:,j),5,2,1);
   % Rotate estimate
   rho     = corr(cumsum([0;du(:,1)]),utcsa(:,j));
   if rho < 0; 
       du  = -du; 
       dl  = -dl; 
   end;
   ufac    = cumsum([zeros(1,size(du,2));du]);
   dufac(:,:,j) = du;
   dlam(:,:,j)  = dl;
   deig(:,j)  = dv;
   % Calibrate to cross-section mean
   sd      = std(utcsa(:,j));
   mn      = mean(utcsa(:,j));
   p0      = [1,0.5];
   opt     = optimset('tolfun',1e-50,'display','off');
   [p,obj] = fminsearch(@(p)calibratef(p,ufac(:,1),sd,mn),p0,opt);
   utpca(:,j) = exp((p(1)*ufac(:,1)+p(2))./2); 
end

% Save results
save aggu dates ut utpca utcsa dufac dlam deig dlogu



