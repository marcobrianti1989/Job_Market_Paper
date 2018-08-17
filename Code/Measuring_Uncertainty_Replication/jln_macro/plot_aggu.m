% -------------------------------------------------------------------------
% Plot aggregate uncertainty estimates
% -------------------------------------------------------------------------

% Load estimates
clear; clc; close all;
load aggu;
T = length(dates);

% Load raw industrial production
[data,txt] = xlsread('jlnrawdata.xlsx',1);
ipg = data(:,6);
ipg = [NaN;log(ipg(2:end)./ipg(1:end-1))];
ipg = tsmovavg(ipg,'s',12,1);
ipg = ipg(end-T+1:end)*100;

% Plot csa estimates
fig = figure(1);
u1  = utcsa(:,1);
u3  = utcsa(:,3);
u12 = utcsa(:,12);
plot(dates,u1,'b','linewidth',1.5); hold on;
plot(dates,u3,'k','linewidth',1);
plot(dates,u12,'-.r','linewidth',1.5);
plot([dates(1),dates(end)],(mean(u1)+1.65*std(u1)).*[1,1],'--b');
plot([dates(1),dates(end)],(mean(u3)+1.65*std(u3)).*[1,1],'--k');
plot([dates(1),dates(end)],(mean(u12)+1.65*std(u12)).*[1,1],'--r');
txt12 = '$\overline{\mathcal{U}}^y_t(12)$';
txt3  = '$\overline{\mathcal{U}}^y_t(3)$';
txt1  = '$\overline{\mathcal{U}}^y_t(1)$';
text(1993,0.95,txt12,'interpreter','latex','color','r');
text(1984,0.90,txt3,'interpreter','latex','color','k');
text(2002,0.59,txt1,'interpreter','latex','color','b');
lab1 = sprintf('h = 1, corr with IP = %0.2f',corr(u1,ipg));
lab2 = sprintf('h = 3, corr with IP = %0.2f',corr(u3,ipg));
lab3 = sprintf('h = 12, corr with IP = %0.2f',corr(u12,ipg));
leg = legend(lab1,lab2,lab3);
set(leg,'location','northwest','box','off');
xlim([dates(1),dates(end)]);
ylim([0.5,1.5]);
rshade(dates);

% Print figure
dim = [6,5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','aggu_csa');

% Plot pca estimates
fig = figure(2);
u1  = utpca(:,1);
u3  = utpca(:,3);
u12 = utpca(:,12);
plot(dates,u1,'b','linewidth',1.5); hold on;
plot(dates,u3,'k','linewidth',1);
plot(dates,u12,'-.r','linewidth',1.5);
plot([dates(1),dates(end)],(mean(u1)+1.65*std(u1)).*[1,1],'--b');
plot([dates(1),dates(end)],(mean(u3)+1.65*std(u3)).*[1,1],'--k');
plot([dates(1),dates(end)],(mean(u12)+1.65*std(u12)).*[1,1],'--r');
txt12 = '$\widehat{\mathcal{U}}^y_t(12)$';
txt3  = '$\widehat{\mathcal{U}}^y_t(3)$';
txt1  = '$\widehat{\mathcal{U}}^y_t(1)$';
text(1993,0.96,txt12,'interpreter','latex','color','r');
text(1984,0.92,txt3,'interpreter','latex','color','k');
text(2002,0.56,txt1,'interpreter','latex','color','b');
lab1 = sprintf('h = 1, corr with IP = %0.2f',corr(u1,ipg));
lab2 = sprintf('h = 3, corr with IP = %0.2f',corr(u3,ipg));
lab3 = sprintf('h = 12, corr with IP = %0.2f',corr(u12,ipg));
leg = legend(lab1,lab2,lab3);
set(leg,'location','northwest','box','off');
xlim([dates(1),dates(end)]);
ylim([0.5,1.5]);
rshade(dates);

% Print figure
dim = [6,5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','aggu_pca');