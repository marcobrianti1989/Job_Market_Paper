clear
close all

figure(1)
set(gcf,'color','w');
plot([0 0 0],[0 1 2],'linewidth',2,'color','k')
xlim([0 1])
ylim([0 2])
yticks([0 1 2])
yticklabels({'Period 2','Period 1','Period 0'})
xticks([])
box off