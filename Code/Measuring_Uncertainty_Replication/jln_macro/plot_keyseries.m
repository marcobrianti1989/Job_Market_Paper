% -------------------------------------------------------------------------
% Compare uncertainty estimates for key series
% -------------------------------------------------------------------------

% Load data
clear; clc;
load ferrors; % for variable names
load npaggu;
h       = 1; % horizon to consider
ut0     = ut(:,:,h);
load aggu;
ut1     = ut(:,:,h);
[T,N,h] = size(ut);
ut0     = ut0(end-T+1:end,:);

% Select series to compare
ind     = [6,37,51,132,72,115,100,95];

% Plot results
fig = figure(1);
for i = 1:length(ind);
   subplot(4,2,i);
   plot(dates,sqrt(ut1(:,ind(i))),'b','linewidth',1.5); hold on;
   plot(dates,sqrt(ut0(:,ind(i))),'k'); hold off;
   title(names(ind(i)));
   xlim([dates(1),dates(end)]);
   if i<6; ylim([0,3]); end;
   if i>6; ylim([0,4]); end;
   if i==3;
       leg = legend('Baseline','No predictors');
       set(leg,'location','northwest','box','off');
   end
   rshade(dates);
end

% Print figure
dim = [6,7];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','keyseries');