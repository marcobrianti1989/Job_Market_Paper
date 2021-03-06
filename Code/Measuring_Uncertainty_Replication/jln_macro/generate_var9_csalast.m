% -------------------------------------------------------------------------
% Estimate an 9-variable 12-th order vector autoregression
% -------------------------------------------------------------------------

% Initialization
clear; clc;
filename = 'var9_csalast';
for ii = [1,2,3,0];
estimate = ii;
bound    = 0; % confidence bound for variance decompositions?
windows=1; mac=0;

% Load and prepare data
if estimate == 1;
load aggu;
[T,N,h]    = size(ut);
[data,txt] = xlsread('vardata.xlsx');
if mac==1;
    data(:,1)=[];
end;
data       = data(end-T+1:end,:);
lipm       = log(data(:,1));
lempm      = log(data(:,2));
hoursm     = data(:,3);
lcpi       = log(data(:,4));
lwage      = log(data(:,5));
ffr        = data(:,6);
lstock     = log(data(:,7));
vxo        = data(:,8);
z          = [lstock,ffr,lwage,lcpi,hoursm,lempm,lipm];
[dum,zf]   = hpfilter(z,129600);

% Estimate vector autoregression
for j = 1:h
    u = utcsa(:,j);
    z = [zf,vxo,u];
    [B,e,SIG,Bs,SIGs] = var_bab(z,12,1000,2000,0,1,1);
    P     = chol(SIG,'lower');
    sigs  = diag(P);
    Om    = diag(diag(P));
    Ainv  = tril(P/Om);
    k     = 12*3+1;
    ci    = 0.6827;
    shock = zeros(1,size(P,1));
    shock(end) = 4.124746507071661*sigs(end);
    [ir,iru,ird] = var_ir(Ainv,B,Bs,SIGs,shock,k,ci);
    IR{j}  = ir.*100;
    IRu{j} = iru.*100;
    IRd{j} = ird.*100;
    params{j}.Ainv = Ainv;
    params{j}.B    = B;
    params{j}.Bs   = Bs;
    params{j}.SIGs   = SIGs;
    params{j}.Om   = Om;
end
save(filename,'params','IR','IRu','IRd','k','ci','h');
end

% Estimate variance decompositions
if estimate == 2;
load(filename);
VD  = cell(h,1);
VDu = cell(h,1);
VDd = cell(h,1);
for j = 1:h
    tic;
    Ainv   = params{j}.Ainv;
    B      = params{j}.B;
    Om     = params{j}.Om;
    if bound == 1; Bs = params{j}.Bs; SIGs = params{j}.SIGs; end;
    if bound == 0; Bs = 0; SIGs = 0; end;
    [vd,vdu,vdd] = var_vd(Ainv,B,Bs,SIGs,Om,k,ci,1);
    VD{j}  = vd.*100;
    VDu{j} = vdu.*100;
    VDd{j} = vdd.*100;
    fprintf('Horizon %d, Elapsed Time = %0.4f \n',j,toc);
end
save(filename,'VD','VDu','VDd','-append');
end

% Estimate variance decomposition at infinite horizon
if estimate == 3;
load(filename);
VDi  = cell(length(h),1);
VDui = cell(length(h),1);
VDdi = cell(length(h),1);
for j = 1:h
    tic;
    Ainv            = params{j}.Ainv;
    B               = params{j}.B;
    Om              = params{j}.Om;
    if bound == 1; Bs = params{j}.Bs; SIGs = params{j}.SIGs; end;
    if bound == 0; Bs = 0; SIGs = 0; end;
    [vdi vdui vddi] = var_vd(Ainv,B,Bs,SIGs,Om,1000,ci,0);
    VDi{j}          = vdi.*100;
    VDui{j}         = vdui.*100;
    VDdi{j}         = vddi.*100;
    fprintf('Horizon %d, Elapsed Time = %0.4f \n',j,toc);
end
save(filename,'VDi','VDui','VDdi','-append');
end

% Plot results
if estimate == 0;
load(filename);
x   = 0:k-1;
ind = [7,6,5];
hor = [1,3,12];
c   = {'b','k','r',[0,0.5,0]};
dg  = [0.80,0.80,0.80];
fig = figure(1);
set(gcf,'defaultlinelinewidth',1.5);
for i = 1:length(hor);
    ir = IR{hor(i)}; iru = IRu{hor(i)}; ird = IRd{hor(i)};
    % Production
    subplot(3,2,(i-1)*2+1);
    plot(x,[iru(:,ind(1)),ird(:,ind(1))],'color',c{i},'linestyle','--');
    hold on; xlim([x(1),x(end)]); ylim([-3,3]);
    plot(xlim,[0,0],'k','linewidth',0.5);
    plot(x,ir(:,ind(1)),'color',c{i}); hold off;
    if i==1; title('Production'); end;
    if hor(i)<=12; ylabel(sprintf('h = %0.f',hor(i))); end;
    if hor(i)>12; ylabel('VXO Index'); end;
    % Employment
    subplot(3,2,(i-1)*2+2);
    plot(x,[iru(:,ind(2)),ird(:,ind(2))],'color',c{i},'linestyle','--');
    hold on; xlim([x(1),x(end)]);ylim([-2,2]);
    plot(xlim,[0,0],'k','linewidth',0.5);
    plot(x,ir(:,ind(2)),'color',c{i}); hold off;
    if i==1; title('Employment'); end;
end

% Print figure
dim = [6,6];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','irf9_csalast');

% Display variance decompositions
num = 9; % which variable to consider
for i = 1:length(ind);
    for j = 1:length(VD);
        vd       = VD{j};
        vdi      = VDi{j};
        col      = [squeeze(vd(ind(i),num,:));vdi(ind(i),num)];
        [m,n]    = max(col);
        v(:,j,i) = [col(1),col(3),col(12),col(end),n,m]';
    end
end
chead = ['      u(1)','  u(3)','  u(12)','  vxo'];
rhead = ['h=1  ';'h=3  ';'h=12 ';'h=inf';'maxh ';'h=max'];
disp('Table Format:');
disp(chead);
disp(rhead); disp(' ');
for i = 1:length(ind);
out = num2str(v(:,:,i),'  %.2f');
disp(['Variable ',num2str(ind(i))]);
disp(out);
disp(' ');
end
end
end