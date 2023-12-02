
clc
clear
close all

step1 = 1; % /10

dir_name = 'figures'; mkdir(dir_name);
figure; hold on; 
FM = [.2 .4 .8 .6];
for i=1:4; 
    FinalMolarity=FM(i); m0=0; Profile_duration = 25; 
    um = (1+FinalMolarity-m0)^(1/Profile_duration);
    t=0:0.0001:200; 
    NaCl = -1 + um.^t; 
    if i==4
        plot(t,  NaCl, '--b', 'LineWidth', 1);
    else
        plot(t,  NaCl, '--','color', [0 0 1 .5], 'LineWidth', .5);
    end
end

plot([0 25],[.6 .6], ':k'); plot([25 25],[0 .6], ':k'); 
indx02M = find(round(NaCl,4)==step1/10); indx02M = indx02M(1); 
INDs = indx02M*[1:5]; 

 plot(t(INDs),NaCl(INDs), 'ro');
 
 plot([0 t(indx02M*[floor(1:.5:5.5)])], [0 0 NaCl(indx02M*[floor(1:.5:5)])], 'k', 'LineWidth', 2); box on; 


%%
stairs = [t(INDs)' NaCl(INDs)'];
hstairs = diff([0; stairs(:,2)]); 

for i=1:5
    text(stairs(i,1)-5, stairs(i,2), num2str(round(stairs(i,2),4)))
    text(3, 1+i*.1, ['\Delta', num2str(i), '   ', num2str(round(hstairs(i),4))]); 
end
ylim([0 1.6]); 
yticks([0:.2:1.6]); 
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'defaultTextFontSize',12, 'FontName', 'Helvetica');

figname = [dir_name, '/stairs_0', num2str(step1), 'M'];  
print(figname,'-depsc', '-r600');
