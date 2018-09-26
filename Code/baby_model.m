clear 
close all

% Calibration
Ushock   = -0.1;
Fshock   = -0.9;
D        = 1;
DU       = D + Ushock;
del1     = 0;
del1F    = del1 + Fshock;
del2     = 0;
d1       = linspace(1.6,2.5,1000);
d2       = linspace(1.6,2.5,1000);

% Model Solve - fmincon
X0       = [1 1];
% SS solve
fun      = @(X) (  D*(1 + (X(2)/X(1)).^0.5) - del1 - X(1)  )^2  ...
      + (   D*(1 + (X(1)/X(2)).^0.5) - del2 - X(2)  )^2  ...
      + (  (  1/(X(1) + del1) + 1/(X(2) + del2)  )^(-1) - D   )^2;  
XX       = fminsearch(fun,X0);
P1opt    = XX(1);
P2opt    = XX(2); 
Q1opt    = P1opt + del1;
Q2opt    = P2opt + del2;

% Financial shock solve
funF     = @(X) (  D*(1 + (X(2)/X(1)).^0.5) - del1F - X(1))^2  ...
      + (   D*(1 + (X(1)/X(2)).^0.5) - del2 - X(2)  )^2  ...
      + (  (  1/(X(1) + del1F) + 1/(X(2) + del2)  )^(-1) - D   )^2;
XXF      = fminsearch(funF,X0);
P1optF   = XXF(1);
P2optF   = XXF(2); 
Q1optF   = P1optF + del1F;
Q2optF   = P2optF + del2;

% Uncertainty shock solve
funU     = @(X) (  DU*(1 + (X(2)/X(1)).^0.5) - del1 - X(1))^2  ...
      + (   DU*(1 + (X(1)/X(2)).^0.5) - del2 - X(2)  )^2  ...
      + (  (  1/(X(1) + del1) + 1/(X(2) + del2)  )^(-1) - DU   )^2;
XXU      = fminsearch(funU,X0);
P1optU   = XXU(1);
P2optU   = XXU(2);
Q1optU   = P1optU + del1;
Q2optU   = P2optU + del2;

% Model Plot in ss
P1d      = P2opt*(d1./D - 1).^(-2);
P2d      = P1opt*(d2./D - 1).^(-2);
P1s      = d1 - del1;
P2s      = d2 - del2;

% Model Plot after F shock
P1dF      = P2optF*(d1./D - 1).^(-2);
P2dF      = P1optF*(d2./D - 1).^(-2);
P1sF      = d1 - del1F;
P2sF      = d2 - del2;

% Model Plot after U shock
P1dU      = P2optU*(d1./DU - 1).^(-2);
P2dU      = P1optU*(d2./DU - 1).^(-2);
P1sU      = d1 - del1;
P2sU      = d2 - del2;

% Figure parameters
sizeline = 3;
gr       = [.7 .7 .7];
lb       = [0.5843 0.8157 0.9882];
% Financial Shock
figure(1)
set(gcf,'Position',[1 41 1920 963])
set(gcf,'color','w');
% External Finance Market
subplot(1,2,1)
plot(d1,P1d,'linewidth',sizeline,'color','k')
hold on
plot(d1,P1s,'linewidth',sizeline,'color',gr)
plot(d1,P1dF,'linewidth',sizeline,'color','b')
plot(d1,P1sF,'linewidth',sizeline,'color',lb)
plot([Q1opt Q1opt],[0 P1opt],'--','color','k')
plot([Q1optF Q1optF],[0 P1optF],'--','color','k')
LGD = legend('demand','supply','demand F Shock','supply F Shock');
LGD.FontSize = 24;
legend boxoff
axis tight
xlabel('Quantity','fontsize',16)
ylabel('Price','fontsize',16)
title('External Sources of Finance','fontsize',24)
set(gca,'YTick',[]);
set(gca,'XTick',[]);
% Internal Finance Market
subplot(1,2,2)
plot(d1,P2d,'linewidth',sizeline,'color','k')
hold on
plot(d1,P2s,'linewidth',sizeline,'color',gr)
plot(d1,P2dF,'linewidth',sizeline,'color','b')
plot([Q2opt Q2opt],[0 P2opt],'--','color','k')
plot([Q2optF Q2optF],[0 P2optF],'--','color','k')
LGD = legend('demand','supply','demand F Shock');
LGD.FontSize = 24;
legend boxoff
axis tight
xlabel('Quantity','fontsize',16)
ylabel('Price','fontsize',16)
title('Internal Sources of Finance','fontsize',24)
set(gca,'YTick',[]);
set(gca,'XTick',[]);
plot([Q2opt Q2optF],[0.2 0.2],'color','r','linewidth',sizeline)
plot([Q2optF-0.05 Q2optF],[0.3 0.2],'color','r','linewidth',sizeline)
plot([Q2optF-0.05 Q2optF],[0.1 0.2],'color','r','linewidth',sizeline)

% Uncertainty Shock
figure(2)
set(gcf,'Position',[-1919 41 1920 963])
set(gcf,'color','w');
% External Finance Market
subplot(1,2,1)
plot(d1,P1d,'linewidth',sizeline,'color','k')
hold on
plot(d1,P1s,'linewidth',sizeline,'color',gr)
plot(d1,P1dU,'linewidth',sizeline,'color','b')
plot([Q1opt Q1opt],[0 P1opt],'--','color','k')
plot([Q1optU Q1optU],[0 P1optU],'--','color','k')
LGD = legend('demand','supply','demand U Shock');
LGD.FontSize = 24;
legend boxoff
axis tight
xlabel('Quantity','fontsize',16)
ylabel('Price','fontsize',16)
title('External Sources of Finance','fontsize',24)
set(gca,'YTick',[]);
set(gca,'XTick',[]);
% Internal Finance Market
subplot(1,2,2)
plot(d1,P2d,'linewidth',sizeline,'color','k')
hold on
plot(d1,P2s,'linewidth',sizeline,'color',gr)
plot(d1,P2dU,'linewidth',sizeline,'color','b')
plot([Q2opt Q2opt],[0 P2opt],'--','color','k')
plot([Q2optU Q2optU],[0 P2optU],'--','color','k')
LGD = legend('demand','supply','demand U Shock');
LGD.FontSize = 24;
legend boxoff
axis tight
xlabel('Quantity','fontsize',16)
ylabel('Price','fontsize',16)
title('Internal Sources of Finance','fontsize',24)
set(gca,'YTick',[]);
set(gca,'XTick',[]);
plot([Q2opt Q2optU],[0.2 0.2],'color','r','linewidth',sizeline)
plot([Q2optU+0.05 Q2optU],[0.3 0.2],'color','r','linewidth',sizeline)
plot([Q2optU+0.05 Q2optU],[0.1 0.2],'color','r','linewidth',sizeline)


% base_path = pwd;
% invoke_export_fig('baby_model_uncertainty_shocks','',0, base_path)


