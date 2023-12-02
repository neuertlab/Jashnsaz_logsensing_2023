path0 = '/Users/hosseinjashnsaz/Dropbox (VU Basic Sciences)/Hossein Jashnsaz/Experiments/Diverse_kinetics_Hog1YFP_TimeLapse2/#7 data structure/';
DATA = load([path0, 'Hog1SignalingData.mat']); 
Hog1SignalingData = DATA.Hog1SignalingData; 
Hog1SignalingData2=Hog1SignalingData; 
%% 
dir_name = 'plots'; mkdir(dir_name);
close all

figure(1); set(gcf, 'Units', 'centimeters', 'Position', [0 0 24 16], 'PaperUnits', 'centimeters', 'PaperSize', [24 16]); 
dx = 0.035; dy = 0.04; 
dd = 6; lw = 2; 

cols = [0 0 0; 0 0 1; 1 0 1]; 

% stimuli 
data = stimuli_functions(); 
subplotHJ(3,4,1,dy,dx); hold on; grid on;
set(groot,'defaultAxesTickLabelInterpreter','tex');
shifts = [0 12 20]; lww=[5 4 2]; lt={'-',':','-','-','-'};  
count = 1; 
% plot(data{1}.tt,zeros(size(data{1}.tt)),'LineWidth', lw, 'color','k');
for i=[6 6 6]    
    plot(data{i}.tt+shifts(count),data{i}.stimuli_profile,'LineWidth', lww(count), 'color',cols(count,:)); 
    count = count + 1; 
end
legend({'3 pre-stimuli yfp','15 pre-stimuli yfp','23 pre-stimuli yfp'}); legend('boxoff')
box on; xlim([-3.5 50.5]); ylim([-.01 .61]); xticks([0:10:50]); yticks([0:.2:.6]); 
ylabel('NaCl (M)'); 

% volume
subplotHJ(3,4,5,dy,dx); hold on; grid on; 
exps = [3 4 5]; dd = 1; 
shift_vol = [0 .03195-.00822 .05661-.00822];  
for i=1:3 
    exp=exps(i); 
    [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
    shadedErrorBar(Hog1SignalingData(exp).tt(dp),-shift_vol(i)+Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cols(i,:)}, 0.2); 
end
box on; xlim([-1.5 53.5]); ylim([.7 1.25]); 
ylabel('Relative volume change');

% Hog1
subplotHJ(3,4,9,dy,dx); hold on; grid on; 
dd = 6;
shifts = [0 0 -1]; 
for i=1:3 
    exp=exps(i);
    [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
	shadedErrorBar(Hog1SignalingData(exp).tt(dp)+shifts(i),Hog1SignalingData(exp).Hog1m(dp),Hog1SignalingData(exp).Hog1s(dp),{'LineWidth', lw, 'color',cols(i,:)}, 0.2); 
end
box on; xlim([-1.5 53.5]); ylim([-.05 0.75]); yticks([0:.1:.7]); 
xlabel('Time (min)'); ylabel('Hog1 nuclear localization');  


%% step 0.6M 7 min with varying pre-stimuli yfp 
exp=6; n_BRs = 9; Hog1SignalingData=Hog1SignalingData2; 

% stimuli
cols = copper(10);
subplotHJ(2,4,2,dy,dx); hold on; grid on;
set(groot,'defaultAxesTickLabelInterpreter','tex');
count = 1; 
% plot(data{1}.tt,zeros(size(data{1}.tt)),'LineWidth', lw, 'color','k');
lww = linspace(6,1,9); 
for i=6*ones(1,9)         
    plot(data{i}.tt,data{i}.stimuli_profile,'LineWidth', lww(count), 'color',cols(count,:)); 
    count = count + 1; 
end
legend({'0 yfp','10 yfp','15 yfp','20 yfp','30 yfp','40 yfp','50 yfp', '60 yfp', '75 yfp'}); 
box on; xlim([-0.5 7.5]); ylim([-.01 .61]); xticks([0:1:7]); yticks([0:.2:.6]); 
ylabel('NaCl (M)'); 

% % volume
% subplotHJ(3,3,6,dy,dx); hold on; grid on; 
% % n_BRs = size(Hog1SignalingData(exp).Sel,2); 
% [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData); 
% dpV=[1:dp(1)-1 dp(1):1:dp(end)]; 
% for br=1:n_BRs
%     plot(Hog1SignalingData(exp).BRs_t(dpV,br), Hog1SignalingData(exp).BRs_VOL(dpV,br),'LineWidth',2, 'color',cols(br,:))
% end
% box on; xlim([-1.2 6.5]); ylim([.7 1.05]);

% Hog1 
[Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData); 
subplotHJ(2,4,6,dy,dx); hold on; grid on; 
dpH=[1:dp(1)-1 dp(1):3:dp(end)]; 
smoothing = false;  stdplot=false;
for br=1:n_BRs
    if stdplot
    if smoothing
        for c=1:size(Hog1SignalingData(exp).Sel(br).InIcSelNintpSBack2N,1)
            cells(c,:) = smooth(Hog1SignalingData(exp).Sel(br).InIcSelNintpSBack2N(c,end-43:end));
        end
    else
        cells = Hog1SignalingData(exp).Sel(br).InIcSelNintpSBack2N(:,end-43:end);
    end
    plotCI((-8:35)/6, cells, cols(br,:)); 
    end
    if smoothing
        plot(Hog1SignalingData(exp).BRs_t(dpH,br), smooth(Hog1SignalingData(exp).BRs_Hog1m(dpH,br)),'-o','LineWidth',2, 'color',cols(br,:))
    else
        plot(Hog1SignalingData(exp).BRs_t(dpH,br), Hog1SignalingData(exp).BRs_Hog1m(dpH,br),'-o','LineWidth',2, 'color',cols(br,:))
    end
end
box on; xlim([-1.2 6.5]); ylim([-.05 0.75]); xticks([0:1:7]); yticks(.1*[0:1:7]); 
xlabel('Time (min)'); ylabel('Hog1 nuclear localization');  

%% Hog1nuc max vs yfp and fit exponential decay
exp=6; 
MaxHog1=max(Hog1SignalingData(exp).BRs_Hog1m,[],1); 
cmap = copper(10);
subplotHJ(2,4,3,dy,dx); hold on; grid on;
set(groot,'defaultAxesTickLabelInterpreter','tex');
count = 1; 
fyp_times=[0 10 15 20:10:60 75];
for i=1
    plot(fyp_times,MaxHog1(i,:),'-','color',[0 0 0]+i*.05); 
    for j=1:9
    plot(fyp_times(j),MaxHog1(i,j),'*','color',cmap(j,:))
    end
    %% exp2 fit
    [MaxHog1_fits(i,:),ff{i}.fit_fun] = expfit(fyp_times,MaxHog1(i,:));
%     plot(fyp_times,MaxHog1_fits(i,:),'-','color',[0 1 1 2/(i*.2+2)]);
    plot([0:.5:100],ff{i}.fit_fun(0:.5:100),'-','color',[0 1 1 2/(i*.2+2)]); 
end
box on; xlim([-0.5 100.5]); ylim([-.01 .71]); xticks([0:5:75 100]); yticks([0:.1:.6]); 
ylabel('maximum(Hog1nuc)');  


%% 10 max Hog1nuc vs yfp and fit exponential decay
Hog1SignalingData(exp).BRs_Hog1m(isnan(Hog1SignalingData(exp).BRs_Hog1m))=0; HOG1nuc=sort(Hog1SignalingData(exp).BRs_Hog1m,'descend'); 
MaxHog1 =  HOG1nuc(1:3:30,:); 

cmap = copper(10);
subplotHJ(2,4,7,dy,dx); hold on; grid on;
set(groot,'defaultAxesTickLabelInterpreter','tex');
count = 1; 
fyp_times=[0 10 15 20:10:60 75];
for i=1:10
    plot(fyp_times,MaxHog1(i,:),'-','color',[0 0 0]+i*.05); 
    for j=1:9
    plot(fyp_times(j),MaxHog1(i,j),'*','color',cmap(j,:))
    end
    %% exp2 fit
    [MaxHog1_fits(i,:),ff{i}.fit_fun] = expfit(fyp_times,MaxHog1(i,:));
%     plot(fyp_times,MaxHog1_fits(i,:),'-','color',[0 1 1 2/(i*.2+2)]);
    plot([0:.5:100],ff{i}.fit_fun(0:.5:100),'-','color',[0 1 1 2/(i*.2+2)]); 
end
box on; xlim([-0.5 100.5]); ylim([-.01 .71]); xticks([0:5:75 100]); yticks([0:.1:.6]); 
xlabel('# of yfp, Time (min)'); ylabel('10 maximum(Hog1nuc)');  

Hog1nuc_scaling.fyp_times = fyp_times;  
Hog1nuc_scaling.MaxHog1 = MaxHog1;  
Hog1nuc_scaling.MaxHog1_fits = MaxHog1_fits;  
Hog1nuc_scaling.MaxHog1_scaling = MaxHog1_fits./(MaxHog1_fits(:,1));  
Hog1nuc_scaling.ff = ff; 
% save('Hog1nuc_scaling','Hog1nuc_scaling'); 

%% correct Hog1 for step 0.6M 7 min with varying pre-stimuli yfp
subplotHJ(3,4,4,dy,dx); hold on; grid on; 
dd=3; 
smoothing = false;  stdplot=true;

[Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData); 

for br=1:n_BRs
    if stdplot
    if smoothing
        for c=1:size(Hog1SignalingData(exp).Sel(br).InIcSelNintpSBack2N,1)
            cells(c,:) = smooth(Hog1SignalingData(exp).Sel(br).InIcSelNintpSBack2N(c,end-43:end));
        end
    else
        cells = Hog1SignalingData(exp).Sel(br).InIcSelNintpSBack2N(:,end-43:end);
    end
    plotCI((-8:35)/6, cells, cmap(br,:)); 
    else
    if smoothing
        plot(Hog1SignalingData(exp).BRs_t(dp,br), smooth(Hog1SignalingData(exp).BRs_Hog1m(dp,br)),'-o','LineWidth',2, 'color',cmap(br,:))
    else
        plot(Hog1SignalingData(exp).BRs_t(dp,br), Hog1SignalingData(exp).BRs_Hog1m(dp,br),'-o','LineWidth',2, 'color',cmap(br,:))
    end
    end
    scalor = [];  
    for i=1:10
        scalor(:,i) = ff{i}.fit_fun(0)./ff{i}.fit_fun(fyp_times(br)+Hog1SignalingData(exp).BRs_t(dp,br));
    end
    mscalor = nanmean(scalor,2); size(mscalor); 
%     mscalor=scalor(1,:); 
    plot(Hog1SignalingData(exp).BRs_t(dp,br), mscalor.*Hog1SignalingData(exp).BRs_Hog1m(dp,br),'-o','LineWidth',2, 'color',[0 0 1 .5]); 
end
box on; xlim([-1.2 6.5]); ylim([-.05 0.85]); xticks([0:1:7]); yticks(.1*[0:1:7]); 
ylabel('Hog1 nuclear localization'); 

%% Hog1nuc and corrected Hog1nuc for 0min, 14min, 20min
path0 = '/Users/hosseinjashnsaz/Dropbox (VU Basic Sciences)/Hossein Jashnsaz/Experiments/Diverse_kinetics_Hog1YFP_TimeLapse2/#7 data structure/';
DATA = load([path0, 'Hog1SignalingData_steps06M.mat']);
Hog1SignalingData = DATA.Hog1SignalingData;
Hog1SignalingData2=Hog1SignalingData; 

load Hog1nuc_scaling
ff = Hog1nuc_scaling.ff; 

lw=2; cols = [0 0 0; 0 0 1; 1 0 1]; 

step_shifts=[0 14 20];  
meths={'meth1','meth2','meth3', 'meth4', 'meth5', 'meth6', 'meth7', 'meth8','meth9','meth10'};
cmap = hsv(length(meths)); 
dd=6; 
scaled_mHog1nuc_tot = []; 
shifts = [0 0 -1]; 
for meth= 8 %1:length(meths)
    
    % measured & analysed data (corrected for single cells photobleaching)
    subplotHJ(3,4,8,dy,dx); hold on; grid on;
    set(groot,'defaultAxesTickLabelInterpreter','tex');
%     title(meths{meth},'color',cmap(meth,:))
    for j=[1:3]
        exp = 3*(meth-1)+(j);
        [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData); 
        shadedErrorBar(Hog1SignalingData(exp).tt(dp)+shifts(j),Hog1SignalingData(exp).Hog1m(dp),Hog1SignalingData(exp).Hog1s(dp),{'LineWidth', lw, 'color',cols(j,:)}, 0.2); 
    end
    box on; xlim([-1.5 50.5]); ylim([-.02 1.52]); xticks([0:10:50]); 
    ylabel('Hog1 nuclear localization');  
    
    % scaled data
    subplotHJ(3,4,12,dy,dx); hold on; grid on;
    set(groot,'defaultAxesTickLabelInterpreter','tex');
    scaled_mHog1nuc = NaN(3,1e3); 
    for j=[1:3]
        exp = 3*(meth-1)+(j);
        [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData); 
        
        % calculate scalor from the 7mins step data 
        scalor = [];  
        for i=1:4
            scalor(i,:) = ff{i}.fit_fun(0)./min(ff{i}.fit_fun(0),ff{i}.fit_fun(step_shifts(j)+Hog1SignalingData(exp).tt(dp)));
        end
        max(scalor); 
        mscalor = nanmean(scalor); size(mscalor); 
        scaled_mHog1nuc(j,1:length(mscalor)) = mscalor.*Hog1SignalingData(exp).Hog1m(dp); 
        shadedErrorBar(Hog1SignalingData(exp).tt(dp)+shifts(j),scaled_mHog1nuc(j,1:length(mscalor)),Hog1SignalingData(exp).Hog1s(dp),{'LineWidth', lw, 'color',cols(j,:)}, 0.2); 
    end
    box on; xlim([-1.5 50.5]); ylim([-.02 1.52]); xticks([0:10:50]); 
    xlabel('Time (min)'); ylabel('Hog1 nuclear localization');  

    scaled_mHog1nuc_tot = [scaled_mHog1nuc_tot; scaled_mHog1nuc(1,:) - scaled_mHog1nuc(2,:) scaled_mHog1nuc(1,:) - scaled_mHog1nuc(3,:) scaled_mHog1nuc(2,:) - scaled_mHog1nuc(3,:)];  
end

%% save figures
set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica')
figname = [dir_name, '/FigureS02_ShiftedSteps_and_Photobleach_Correction'];  
print(figname,'-depsc', '-r600');


