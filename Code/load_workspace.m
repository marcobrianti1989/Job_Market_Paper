clear
close all

load workspace_chol_Ufirst_GZSpread
vardec_Ufirst = vardec; 
ss_Ufirst     = ss(:,1);

% plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,which_shocks,shocknames,...
%       system_names,which_ID,print_figs,use_current_time,base_path)

IRFsA11 = IRFs(:,:,1);
IRFsA111 = [IRFsA11(2:end,:); IRFsA11(1,:)];
IRFsA = zeros(size(IRFs,1),size(IRFs,2),size(IRFs,3));
IRFsA(:,:,end) = IRFsA111;

ubA11 = ub2(:,:,1);
ubA111 = [ubA11(2:end,:); ubA11(1,:)];
ubA = zeros(size(IRFs,1),size(IRFs,2),size(IRFs,3));
ubA(:,:,end) = ubA111;

lbA11 = lb2(:,:,1);
lbA111 = [lbA11(2:end,:); lbA11(1,:)];
lbA = zeros(size(IRFs,1),size(IRFs,2),size(IRFs,3));
lbA(:,:,end) = lbA111;

load workspace_chol_Ulast_GZSpread
vardec_Ulast  = vardec;
ss_Ulast      = ss(:,end);

plot_single_IRFs_2CIs_2specifications(IRFsA,ubA,lbA,IRFs,ub2,lb2,H,...
      which_shocks,shocknames,...
      system_names,which_ID,print_figs,use_current_time,base_path)

close all

return

%Define Forecast Horizon
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
title('Variance Explained Of Real GDP')

figure(2)
plot(horz,vardec_Ufirst(6,:),'linewidth',2)
hold on
plot(horz,vardec_Ulast(5,:),'linewidth',2)
lgd = legend('Unconstrained Uncertainty Shock','Constrained Uncertainty Shock');
lgd.FontSize = 24;
lgd.Location = 'northwest';
grid on
legend boxoff
xlabel('Horizon')
ylabel('Variance Explained')
title('Variance Explained Of Real Investment')

figure(3)
plot(horz,vardec_Ufirst(5,:),'linewidth',2)
hold on
plot(horz,vardec_Ulast(4,:),'linewidth',2)
lgd = legend('Unconstrained Uncertainty Shock','Constrained Uncertainty Shock');
lgd.FontSize = 24;
lgd.Location = 'northwest';
grid on
legend boxoff
xlabel('Horizon')
ylabel('Variance Explained')
title('Variance Explained Of Real Consumption')

figure(4)
plot(Time(nlags+1:end),ss_Ufirst,'linewidth',2)
hold on
plot(Time(nlags+1:end),ss_Ulast,'linewidth',2)
lgd = legend('Unconstrained Uncertainty Shock','Constrained Uncertainty Shock');
lgd.FontSize = 24;
lgd.Location = 'northwest';
grid on
legend boxoff

corr(ss_Ufirst,ss_Ulast)