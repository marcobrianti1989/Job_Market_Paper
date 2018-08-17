% -------------------------------------------------------------------------
% Estimate 8 variable 12th order vector autoregression for other proxies of
% uncertainty: sd firm stock returns, sd firm profits, sd GDP forecasts,
% and sd industry-level total factor productivity.
% -------------------------------------------------------------------------

% Initialization
clear; clc; close all;
filename = 'cee05var_proxies';
for ii = [0]
estimate = ii;

windows=0; mac=1;
% Load data
if estimate == 1;
load proxies;
[data,txt] = xlsread('ceevardata.xlsx',2);
if mac==1;
    data(:,1)=[];
end;
lip        = log(data(:,1));
lemp       = log(data(:,2));
lcons      = log(data(:,3));
lpr        = log(data(:,4));
lorders    = log(data(:,5))-lpr;
lwage      = log(data(:,6))-lpr;
hrs        = data(:,7);
lhrs       = log(data(:,7));
ffr        = data(:,8)./100;
lstock     = log(data(:,9));
lmzm       = log(data(:,10));
gmzm       = [0;lmzm(2:end)-lmzm(1:end-1)].*1200;
vxo        = data(:,11);
lm1        = log(data(:,12));
gm1        = [0;lm1(2:end)-lm1(1:end-1)].*1200;
lm2        = log(data(:,13));
gm2        = [0;lm2(2:end)-lm2(1:end-1)].*1200;
lcpi       = log(data(:,14));
lcom       = log(data(:,15));
T          = length(data);
dates      = 1900+(50:1/12:112-1/12)';
dates      = dates(end-T+1:end);
xxx
% Stock returns
ind      = find(ret(:,1)==dates(1));
u        = ret(ind:end,2);
z        = [lip,lemp,lcons,lpr,lorders,lwage,lhrs,ffr,lstock,gm2];
zf       = z;
z        = [zf,u];
[B,e,SIG,Bs,SIGs] = var_bab(z,12,1000,5000,0,1,1);
P     = chol(SIG,'lower');
sigs  = diag(P);
Om    = diag(diag(P));
Ainv  = tril(P/Om);
k(1)  = 12*5+1;
ci    = 0.6827;
shock = zeros(1,size(P,1));
shock(end) = 4*sigs(end);
[ir,iru,ird] = var_ir(Ainv,B,Bs,SIGs,shock,k(1),ci);
IR{1}  = ir.*100;
IRu{1} = iru.*100;
IRd{1} = ird.*100;

% Profits
num = isnan(prof(:,2))<1;
yq  = prof(num,1);
u   = prof(num,2);
a   = find(dates==yq(1));
b   = find(dates==yq(end));
z   = zf(a:b+2,:);
z   = tsmovavg(z,'s',3,1);
z   = z(3:3:end,:);
z   = [z,u];
[B,e,SIG,Bs,SIGs] = var_bab(z,4,1000,5000,0,1,1);
P     = chol(SIG,'lower');
sigs  = diag(P);
Om    = diag(diag(P));
Ainv  = tril(P/Om);
k(2)  = 4*5+1;
ci    = 0.6827;
shock = zeros(1,size(P,1));
shock(end) = 4*sigs(end);
[ir,iru,ird] = var_ir(Ainv,B,Bs,SIGs,shock,k(2),ci);
IR{2}  = ir.*100;
IRu{2} = iru.*100;
IRd{2} = ird.*100;

% Forecasts
num = isnan(gdpf2(:,5))<1;
yh  = gdpf2(num,1);
u   = gdpf2(num,5);
a   = find(yh==dates(1));
b   = find(yh<=dates(end),1,'last');
u   = detrend(u(a:b));
z   = zf;
z   = tsmovavg(z,'s',6,1);
z   = z(6:6:end,:);
z   = [z,u];
[B,e,SIG,Bs,SIGs] = var_bab(z,2,1000,5000,0,1,1);
P     = chol(SIG,'lower');
sigs  = diag(P);
Om    = diag(diag(P));
Ainv  = tril(P/Om);
k(3)  = 2*5+1;
ci    = 0.6827;
shock = zeros(1,size(P,1));
shock(end) = 4*sigs(end);
[ir,iru,ird] = var_ir(Ainv,B,Bs,SIGs,shock,k(3),ci);
IR{3}  = ir.*100;
IRu{3} = iru.*100;
IRd{3} = ird.*100;

% Industry productivity
yy = newtfp(:,1);
a  = find(yy>=dates(1),1,'first');
b  = find(dates<=yy(end),1,'last');
u  = newtfp(a:end,2);
z  = zf(1:b+12,:);
z  = tsmovavg(z,'s',12,1);
z  = z(12:12:end,:);
z  = [z,u];
[B,e,SIG,Bs,SIGs] = var_bab(z,1,1000,5000,0,1,1);
P     = chol(SIG,'lower');
sigs  = diag(P);
Om    = diag(diag(P));
Ainv  = tril(P/Om);
k(4)  = 2*5+1;
ci    = 0.6827;
shock = zeros(1,size(P,1));
shock(end) = 4*sigs(end);
[ir,iru,ird] = var_ir(Ainv,B,Bs,SIGs,shock,k(4),ci);
IR{4}  = ir.*100;
IRu{4} = iru.*100;
IRd{4} = ird.*100;

% Save results
save(filename,'IR','IRu','IRd','k');
end

% Plot results
if estimate == 0;
load(filename);
fig = figure(1);
set(gcf,'defaultlinelinewidth',1.5);
i1   = 1;
i2   = 2;
ylab = {'Returns','Profits','Forecasts','TFP'};
xlab = {'Months','Quarters','Half-years','Years'};
for i = 1:length(k)
    x   = 0:k(i)-1;
    subplot(4,2,(i-1)*2+1);
    plot(x,IR{i}(:,i1)); hold on;
    plot(x,IRu{i}(:,i1),'--');
    plot(x,IRd{i}(:,i1),'--');
    xlim([x(1),x(end)]); ylim([-7,7]);
    if i == 1; ylim([-5,5]); end;
    if i == 4; ylim([-10,5]);end;
    plot(xlim,[0,0],'k','linewidth',0.5);
    ylabel(ylab{i}); xlabel(xlab{i});
    if i==1; title('Production'); end;

    subplot(4,2,(i-1)*2+2);
    plot(x,IR{i}(:,i2)); hold on;
    plot(x,IRu{i}(:,i2),'--');
    plot(x,IRd{i}(:,i2),'--');
    xlim([x(1),x(end)]); ylim([-5,5]);
    if i == 4; ylim([-10,5]); end;
    plot(xlim,[0,0],'k','linewidth',0.5);
    if i==1; title('Employment'); end;
end

% Print figure
dim = [6,6];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','cee05irf_proxies');
end
end