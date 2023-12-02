
% load('../../Figure02/analysis and plots meth8 photobleach corr/Figure02'); 

dir_name = 'plots'; mkdir(dir_name);
close all
addpath('../../data/matlab-general-functions/')

%% figure handle
close; figure(4); set(gcf, 'Units', 'centimeters', 'Position', [0 0 5.5 10.5], 'PaperUnits', 'centimeters', 'PaperSize', [5.5 10.5]);
dx = 0.015; dy = 0.03; cmap = jet(7); lw=2; cmap([3 5],:)=.95*cmap([3 5],:); cmap = [cmap; hex2rgb('800080')]; cmap = [0,0,0;cmap]; 
  
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')


%% plot polynomials for stimuli, volume, and Hog1 changes
time=-2.5:.1:27.5;
parameters = [1 NaN 25 0 0]; % (y0,n,T0,T1,T2)
polynomial_orders = [.01 .5 1 2 3 5 7 1]; % t0, t_sqr, t1, t2, t3, t5, t7, exp(t)

load poly_order_vol % volume k fits
load poly_order_Hog1 % Hog1 k fits
for conc = 1:5
for j=1:4
    id = 3+j;
    hh=subplotHJ(4,1,j,dy,dx); cla(hh); hold on; grid on; box on; axis off
    xline([0 25], ':k');
    yline([0 1], ':k');
    
    % stimuli 
    parameters(2) = polynomial_orders(id); 
    stimuli = polynomial_fxn(time,parameters);
%     plot(time,stimuli,'color',.5*[1 1 1], 'LineWidth',2)

     % volume 
    parameters(2) = poly_order_vol(conc,id); 
    volume = polynomial_fxn(time,parameters);
    plot(time,volume,':','color',cmap(id+1,:), 'LineWidth',2)

    % volume 
    parameters(2) = poly_order_Hog1(conc,id); 
    hog1 = polynomial_fxn(time,parameters);
    plot(time,hog1,'color',cmap(id+1,:), 'LineWidth',2)

    xlim([-3 28]); xticks([0:25:50]); xticklabels([]); ylim([0 1]); yticks([0:1:1]); 
end
%% save figs
% set(findall(gcf,'-property','FontSize'),'FontSize',7, 'color',.95*[1 1 1], 'defaultTextFontSize',7, {'DefaultAxesXColor','DefaultAxesYColor'},{.5*[1 1 1],.5*[1 1 1]}, 'FontName', 'Helvetica')
set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica')

figname = [dir_name, '/Figure04_G_',num2str(conc)];
print(figname,'-depsc', '-r600');
end