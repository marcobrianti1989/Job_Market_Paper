% -------------------------------------------------------------------------
% Plot different measures of uncertainty in the S&P 500 index
% -------------------------------------------------------------------------

% Load data
clear; clc;close all;
ind = 82; % S&P 500 index
h   = 1; % horizon to consider
load ut; 
u1     = sqrt(squeeze(ut(:,ind,h)));
T1     = length(u1);
dates1 = dates;
load arut;
u2     = sqrt(squeeze(ut(:,ind,h)));
T2     = length(u2);
dates2 = dates;
load nput;
u3     = sqrt(squeeze(ut(:,ind,h)));
T3     = length(u3);
dates3 = dates;

% Plot data
fig = figure(1);
set(gcf,'defaultlinelinewidth',1.5);
plot(dates1,u1,'b'); hold on;
plot(dates2,u2,'k','linewidth',1);
plot(dates3,u3,'-.r'); hold off;
text(1998.5,2,'No Predictors','color','r');
text(1993,1.4,'AR only','color','k');
annotation('arrow',[0.65,0.61],[0.48,0.3],'color','k','linewidth',1);
text(1975,0.5,'Baseline','color','b');
xlim([dates3(1),dates3(end)]);
rshade(dates3);

% Print figure
dim = [6,5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','usp500');