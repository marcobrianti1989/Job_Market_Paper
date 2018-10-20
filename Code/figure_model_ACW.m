clear
close all

Xlim = 5;
X = linspace(0,Xlim,Xlim*10);

par1 = 1.4;
par2 = 1;
Y = -1 + par1./(par2 + X);

figure(1)
set(gcf,'color','w');
hold on
plot(X,Y,'linewidth',4,'color','k')
plot([0 Xlim],[0 0],'linewidth',1,'color','r')
ylim([-1 1])
xlim([0 Xlim])
set(gca,'XTick',[0],'XTickLabel',{'0'},'fontsize',32);
set(gca,'YTick',[-1 0 +1],'YTickLabel',{'-1','0','1'},'fontsize',32);
xlabel('\delta: Instrument Weight','fontsize',45)
ylabel('\gamma_U \gamma_F'': Shocks Correlation','fontsize',45)
hold off

% base_path         = pwd;
% invoke_export_fig('GPFA_intuition','figure',0,base_path)























asd

figure(1)
set(gcf,'color','w');
plot([0 0 0],[0 1 2],'linewidth',3,'color','k')
xlim([0 1])
ylim([0 2])
set(gca,'YTick',[0 1 2],'YTickLabel',{'Period 2','Period 1','Period 0'},'fontsize',24);
xticks([])
set(gca,'xcolor',[1 1 1]); 
txt0 = 'd_0 = Y_0 + B_0 - I_0 - C';
text(0.1,1.99,txt0,'fontsize',24)
txt1 = 'd_1 = Y_1 + B_1 - I_1 + C';
text(0.1,0.99,txt1,'fontsize',24)
txt2 = 'd_2 = g(I_0) - B_0 + h(I_1) - B_1';
text(0.1,-0.01,txt2,'fontsize',24)
box off