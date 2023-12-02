
dir_name = 'plots'; mkdir(dir_name);
close all
addpath('../../../data/matlab-general-functions/')

%% figure handle
close; figure(1); set(gcf, 'Units', 'centimeters', 'Position', [0 0 10 6], 'PaperUnits', 'centimeters', 'PaperSize', [10 6]); 
dx = 0.05; dy = 0.08; cmap = jet(7); lw=2;

% subplotHJ(1,1,1,dy,dx); 
hold on; grid on;
set(groot,'defaultAxesTickLabelInterpreter','tex');

c_max = 1;  T = 38;

xline([0 T], ':k');
yline([0 c_max], ':k');

syms x
fplot(heaviside(x), [-10, 50],'b', 'LineWidth',7)

stim = @(t,k) (t<=T).*heaviside(t).*(c_max*(t/T).^k) + (t>T).*c_max; 
tt=[-10:1:50]; 
plot(tt, stim(tt,.4), ':k', 'LineWidth',3)
plot(tt, stim(tt,1), ':','color',.5*[1 1 1], 'LineWidth',4)
plot(tt, stim(tt,3), ':k', 'LineWidth',4)

xlim([-10 52]); ylim([-.0 1.01]); xticks([]); yticks([]); 

%% save figs
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

figname = [dir_name, '/stimuli'];  
print(figname,'-depsc', '-r600');
