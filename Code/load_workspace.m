clear
close all

load workspace_chol_Ufirst
vardec_Ufirst = vardec; 
ss_Ufirst     = ss(:,1);

load workspace_chol_Ulast
vardec_Ulast  = vardec;
ss_Ulast      = ss(:,end);

figure
plot(vardec_Ufirst(4,:))
hold on
plot(vardec_Ulast(3,:))

figure
plot(ss_Ufirst)
hold on
plot(ss_Ulast)

corr(ss_Ufirst,ss_Ulast)