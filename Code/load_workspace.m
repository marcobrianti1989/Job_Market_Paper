clear
close all

load workspace_chol_Ufirst_GZSpread
vardec_Ufirst = vardec; 
ss_Ufirst     = ss(:,1);

load workspace_chol_Ulast_GZSpread
vardec_Ulast  = vardec;
ss_Ulast      = ss(:,end);

horz = linspace(0,20,20);

figure(1)
plot(horz,vardec_Ufirst(4,:),'linewidth',2)
hold on
plot(horz,vardec_Ulast(3,:),'linewidth',2)
lgd = legend('Unconstrained Uncertainty Shock','Constrained Uncertainty Shock');
lgd.FontSize = 24;
lgd.Location = 'northwest';
grid on
legend boxoff
xlabel('Horizon')
ylabel('Variance Explained')

figure(2)
plot(ss_Ufirst,'linewidth',2)
hold on
plot(ss_Ulast,'linewidth',2)
lgd = legend('Unconstrained Uncertainty Shock','Constrained Uncertainty Shock');
lgd.FontSize = 24;
lgd.Location = 'northwest';
grid on
legend boxoff

corr(ss_Ufirst,ss_Ulast)