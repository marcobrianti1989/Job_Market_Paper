% -------------------------------------------------------------------------
% Plot uncertainty of predictor variables (h = 1 only)
% -------------------------------------------------------------------------

% Load data
clear; clc; 
load ferrors;
svf = load('svfmeans.txt');

% Compute uncertainty
thf   = [svf(1,:).*(1-svf(2,:));svf(2,:);svf(3,:).^2];
xf    = svf(4:end-3,:);
[T,r] = size(xf);
uf    = zeros(T,r);
for i = 1:r
    uf(:,i) = exp(thf(1,i)+thf(3,i)/2 + thf(2,i).*xf(:,i));
end

% Plot results
fig = figure(1);
set(gcf,'defaultlinelinewidth',1.5);
subplot(3,2,1); plot(dates,uf(:,1)); 
xlim([dates(1),dates(end)]);ylim([0,1.5]);
rshade(dates);
title('$F_{1t}:$ Stock Returns','interpreter','latex');
subplot(3,2,2); plot(dates,uf(:,2)); 
xlim([dates(1),dates(end)]);ylim([0,0.25]);
rshade(dates);
title('$F_{2t}:$ Real Activity','interpreter','latex');
subplot(3,2,3); plot(dates,uf(:,4)); 
xlim([dates(1),dates(end)]);ylim([0,0.25]);
rshade(dates);
title('$F_{4t}:$ Inflation','interpreter','latex');
subplot(3,2,4); plot(dates,uf(:,5)); 
xlim([dates(1),dates(end)]);ylim([0,0.25]);
rshade(dates);
title('$F_{5t}:$ FF Factors and Bond Spreads','interpreter','latex');
subplot(3,2,5); plot(dates,uf(:,13)); 
xlim([dates(1),dates(end)]);ylim([0,4]);
rshade(dates);
title('$F_{1t}^2$','interpreter','latex');
subplot(3,2,6); plot(dates,uf(:,14)); 
xlim([dates(1),dates(end)]);ylim([0,1.5]);
rshade(dates);
title('$G_{1t}$','interpreter','latex');

% Print figure
dim = [6,5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','utpredictor');