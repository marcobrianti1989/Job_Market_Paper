function plot_vardec(IRFs,h,which_shock,names,...
      varnames,which_ID_strat,print_figs,use_current_time,base_path)

% h = IR horizon for figures
% ub1 and lb1 are the bootstrapped CI for the 1st significance level
% ub2 and lb2 are the bootstrapped CI for the 2nd significance level
% names = cell vector of shock names
% varnames = cell vector of variable names
% which_ID_strat = a string describing which identification strategy we used
% print_figs = 'yes' --> saves the figures; else if 'no' -- > just
% shows the figures.

nvar = size(IRFs,1);
nshocks = size(which_shock,2);
periods = 1:h;

% Draw pretty pictures
for i_var=1:nvar
      for i_shock=1:nshocks
            figure(i_var)
            set(gcf,'color','w'); % sets white background color
            set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen (hopefully)
            varname = varnames{i_var};
            name    = names{i_shock};
            hold on
            plot(periods,IRFs(i_var,1:h,which_shock(i_shock)),'linewidth',3)
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', 28)
            title(['Variance Explained of ' , varname],'fontsize',72)
            ylim([0 1]);
            xlim([1 h])
            legendInfo{i_shock} = [names{i_shock}]; % or whatever is appropriate
            LGD = legend(legendInfo);
            LGD.FontSize = 24;
            legend boxoff
            hold off
            grid on
      end
      
      % Save figures if you want to
      if strcmp(print_figs, 'yes')
            invoke_export_fig([name, ' on ' , varname], which_ID_strat,...
                  use_current_time, base_path)
            close all
            pause(0.5)
      end
end














end