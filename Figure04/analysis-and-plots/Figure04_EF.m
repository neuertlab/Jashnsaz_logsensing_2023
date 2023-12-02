%% figure handle
dir_name = 'plots'; mkdir(dir_name);
close all
addpath('../../data/matlab-general-functions/')

close; figure(3); set(gcf, 'Units', 'centimeters', 'Position', [0 0 9 3.5], 'PaperUnits', 'centimeters', 'PaperSize', [9 3.5]);
dx = 0.1; dy = 0.1; cmap = jet(7); lw=2; cmap([3 5],:)=.95*cmap([3 5],:); cmap = [cmap; hex2rgb('800080')]; cmap = [0,0,0;cmap]; 
  
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

%% plot total volume shrink over stimulation time
hh=subplotHJ(1,2,1,dy,dx); cla(hh); hold on; grid on; box on;
set(groot,'defaultAxesTickLabelInterpreter','tex');
        
load total_Vol_shrink % volume k fits
VSh = total_Vol_shrink./total_Vol_shrink(:,2); 
conc=4; 
for i=conc %1:size(VS,1)
    VS=VSh(i,:) 
    for j=1:size(VS,2)
        if j>=2           
            bar(j, VS(j), 'FaceColor',cmap(j,:),'EdgeColor','none');  
        end
    end
end
% VSm=nanmean(VSh); VSs=nanstd(VSh);
% errorbar([1:9],VSm, .5*VSs,.5*VSs,'linewidth',1, 'color', 0*[1 1 1])

xlim([.5 size(VS,2)+.5]); xticks([0:1:size(VS,2)]);xticklabels([]); %ylim([.0 1.01]); yticks([0:.2:1]); % 0:1:N1+1
set(gca,'YScale','log'); ylim([5e-3 1e0]); % yticks([0:.2:1]);

%% Growth rate vs volume shrink
hh=subplotHJ(1,2,2,dy,dx); cla(hh); hold on; grid on; box on;
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(gca,'GridLineStyle',':'); 

% volume reducttion data
VolReduction = total_Vol_shrink./(total_Vol_shrink(:,2)); % each normalized to the step condition
VolReduction_CSLE = VolReduction(:,[1 2 4 end]); % control, step, linear, exponential; all final concentrations 

% cell growth rate data
load GrowthRate_DATA_FINAL_3tp.mat % growth rates from Figure 2 and Figure S1
GrowthRate_SLE = [DoublingRates([1 8 9],:,[2]) DoublingRates([1 8 9],:,[3])]; % pulse 120min, linear, exponential; 3BRs; post-1st-stress and post-2nd-stress
GrowthRate_C = controls; 

errorbar(nanmean(VolReduction_CSLE(:,1)),nanmean(GrowthRate_C), .5*nanstd(GrowthRate_C),.5*nanstd(GrowthRate_C), .5*nanstd(VolReduction_CSLE(:,1)),.5*nanstd(VolReduction_CSLE(:,1)),'linewidth',1, 'color', 0*[1 1 1])
cid = [2 4 9];
for i=1:3
    errorbar(nanmean(VolReduction_CSLE(:,1+i)),nanmean(GrowthRate_SLE(i,:)), .5*nanstd(GrowthRate_SLE(i,:)),.5*nanstd(GrowthRate_SLE(i,:)), .5*nanstd(VolReduction_CSLE(:,1+i)),.5*nanstd(VolReduction_CSLE(:,1+i)),'linewidth',1, 'color', cmap(cid(i),:))
end
xlim([-.05 1.05]); xticks([0:.2:1]); 
ylim([-.05 1.05]); yticks([0:.2:1]);
% set(gca,'XScale','log'); %xlim([5e-4 2e0]); % yticks([0:.2:1]);
% set(gca,'YScale','log'); ylim([5e-2 2e0]); % yticks([0:.2:1]);

%% save fig

set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica')

figname = [dir_name, '/Figure04_EF_',num2str(conc)];
print(figname,'-depsc', '-r600');