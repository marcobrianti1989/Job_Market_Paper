% -------------------------------------------------------------------------
% Compare time series behavior of vxo index and uncertainty estimates
% -------------------------------------------------------------------------

% Load data
clear; clc;
load aggu;
T = length(dates);
%% default is windows
[data,txt] = xlsread('vardata.xlsx');
mac=1; windows=0;
if windows==1;
vxo = data(end-T+1:end,8);
end;
if mac==1;
data(:,1)=[];
txt(:,1)=[];
vxo = data(end-T+1:end,8);
end;


% Dates in Bloom (2009, EMA) Table A.1
bdates = [1962+(10-1)/12,1963+(11-1)/12,1966+(8-1)/12,1970+(5-1)/12,...
          1973+(12-1)/12,1974+(10-1)/12,1978+(11-1)/12,1980+(3-1)/12,...
          1982+(10-1)/12,1987+(11-1)/12,1990+(10-1)/12,1997+(11-1)/12,...
          1998+(9-1)/12,2001+(11-1)/12,2002+(11-1)/12,2003+(2-1)/12,...
          2008+(10-1)/12];

% Threshold computed using updated data
cut     = 1.65;
lam     = 129600;
[tr,cy] = hpfilter(vxo,lam);

% Plot results
gr  = [0,0.5,0];
u   = utcsa(:,1);
fig = figure(1);
plot([dates(1),dates(end)],[cut,cut],'-.r');hold on;
for i = 1:length(bdates);
   plot(bdates(i).*[1,1],[-2,7],'color',gr,'linewidth',0.5);
end
plot(dates,zscore(u),'linewidth',1.5);
plot(dates,zscore(vxo),'k');
txt = '$\overline{\mathcal{U}}^y_t(1)$';
text(1981,4,txt,'interpreter','latex','color','b');
text(1991.5,1.2,'VXO','interpreter','latex','color','k');
text(1962,-1.6,sprintf('Correlation = %0.2f',corr(vxo,u)));
text(1964,2,'1.65 std','color','r');
text(1962,6.5,'Bloom dates','color',gr);
annotation('arrow',[0.29,0.4],[0.875,0.875],'color',gr,'linewidth',0.5);
xlim([dates(1),dates(end)]);
rshade(dates);
hold off;

% Print figure
dim = [6,5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','tsvxo_csa');
