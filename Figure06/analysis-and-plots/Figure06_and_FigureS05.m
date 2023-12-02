
%% load data
warning('all','off')
clear all
%% global data path
% % Hog1SignalingData
% % path0 = '/Users/hosseinjashnsaz/NeuertLab/Dropbox (VU Basic Sciences)/Hossein Jashnsaz/Experiments/Diverse_kinetics_Hog1YFP_TimeLapse/#7 data structure 2 meth8 scale photobleach/'; 
% path0 = '/Users/hosseinjashnsaz/Dropbox (VU Basic Sciences)/Hossein Jashnsaz/Experiments/Diverse_kinetics_Hog1YFP_TimeLapse/#7 data structure 2 meth8 scale photobleach/'; 
% DATA = load([path0, 'Hog1SignalingData.mat']);
% Hog1SignalingData = DATA.Hog1SignalingData;
%  
% % % HogDataSet
% % path0= 'C:\Users\jashnsh\Dropbox (VU Basic Sciences)\Signaling pathway ID\HOGSignalingModeling\ModelData\HJ_20230701_HOGDATA_simple\HogDataSet\';
% % path0= '/Users/hosseinjashnsaz/NeuertLab/Dropbox (VU Basic Sciences)/Signaling pathway ID/HOGSignalingModeling/ModelData/HJ_20230701_HOGDATA_simple/HogDataSet/';
% % DATA = load([path0, 'HogDataSet.mat']);
% % HogDataSet = DATA.HogDataSet;
%% local data path
addpath('../../data/cell-volume-signaling-data/')
load Hog1SignalingData2
load HogDataSet 
dir_name = 'plots'; mkdir(dir_name);
close all
addpath('../../data/matlab-general-functions/')

%% define data
Hog1SignalingData=Hog1SignalingData2; 
DataSet1 = HogDataSet.DataSet1;
DataSet3 = HogDataSet.DataSet3;

% clearvars -except Hog1SignalingData; 
%% set dir
col_csle(1,:) = .6*[1 1 1]; %hex2rgb('f0f5f9'); %ctrl
col_csle(2,:) = hex2rgb('87314e'); %step  090707
col_csle(3,:) = hex2rgb('00b9f1'); %linear
col_csle(4,:) = hex2rgb('A593E0'); %linear
dc = .6; % .27
for ii=1:4
    col = col_csle(ii,:);
    for i=1:4
        cmap{ii}(i,:)=col+(1-col)*dc*(i-1);
    end
end

dd = 6; lw = 1.5; 

Experiments = {'steps $ (t^{0}) $', 'root $ (\sqrt{t}) $', 'linears $ (t^{1}) $', 'quadratics $ (t^{2}) $', 'cubics $ (t^{3}) $', 'quintics $ (t^{5}) $', 'heptics $ (t^{7}) $', 'powers $ (a^t) $', 'sum polys'}; 
ExperimentsID = {[1:9], [10:14], [15:19], [20:24], [25:29], [30:34], [35:39], [40:44]}; 
% ExperimentsID = {[1:3 5 7], [9:13], [14:18], [19:23], [24:28], [29:33], [34:38], [39:43]}; 
finalConcs = {[0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 1.0], [0.10 0.20:0.2:0.8], [[0.10 0.20:0.2:0.8]], ...
    [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]]}; 
%% stimuli 
close; figure();clf; set(gcf, 'Units', 'centimeters', 'Position', [0 0 13.6 9], 'PaperUnits', 'centimeters', 'PaperSize', [13.6 9]); 
dx = 0.008; dy = 0.015; lw=2; 

% load([dir_name,'/stimuli.mat']);
IDS = [29 31; 95 97];

ii=1; 
for i=[1:2]
    subplotHJ(3,3,ii,dy,dx); hold on; %grid on
    count = 1;     
    for exp = 1:2
        id = IDS(i,exp); 
        data = DataSet1{id}; 
        plot(data.tt,-.1+data.stimulus(data.tt),'LineWidth', lw, 'color',cmap{ii+2}(count,:)); 
        box on; xlim([-2.5 32.5]); ylim([-.01 .81]); xticks([0:25:50]); yticks([0:.2:.8]); 
        xticklabels([]); if ii>1; yticklabels([]); end
        count = count + 1; 
    end
    box on; xlim([-2.5 32.5]); ylim([-.01 .81]); xticks([0:25:50]); yticks([0 .6 .8]); 
    xticklabels([]); yticklabels([]); 

    ii=ii+1;
end    

% sum polynomials stimuli

FM = [.2 .4 .6 .8]; % final molarities [.2 .4 .6 .8];
IDs = [3:2:9]; 
polk = [3 4 5 6 7]; 
data = DataSet1{1}; 
sumStimuli = NaN(length(data.stimulus_profile), 4); 

for conc=1:4
    a(conc) = (1+FM(conc)-0)^(1/25);
    count = 1; 
    for polynomials = [1 2 3 5 7] 
        exp = 11*(polk(count)-1)+IDs(conc);
        data = DataSet1{exp}; 
        IN = NaN(length(data.stimulus_profile), 1);
        IN(1:length(data.stimulus_profile),1) = data.stimulus(data.tt);
        sumStimuli(:,conc) = nansum([sumStimuli(:,conc) ((25*log(a(conc))).^polynomials)*IN./(FM(conc)*factorial(polynomials))],2);        
        count  = count + 1; 
    end    
end

subplotHJ(3,3,3,dy,dx); hold on; %grid on
data = DataSet1{1}; 

ii=1;
for i=3:4
    plot(data.tt,-.1+sumStimuli(:,i),'LineWidth', lw, 'color',cmap{2}(ii,:)); 
    ii=ii+1; 
end
box on; xlim([-2.5 32.5]); ylim([-.01 .81]); xticks([0:25:50]); yticks([0 .6 .8]); 
xticklabels([]); yticklabels([]); 


%% volume linear and exponential stimuli
IDS = [18 19; 43 44]; % t1 0.6M 0.8M; exp 0.6M 0.8M

% control 0M
data48 = Hog1SignalingData(48); data48.tt(1)
[time,Im0,Is0] = get_ImIs(data48,-1,25+4,'Vol'); 
size(Im0)
ii=1; 
for i=[1 2]
    hh=subplotHJ(3,3,3+ii,dy,dx); cla(hh); hold on; %grid on
    count = 1; 
    for conc = [1 2]
        data = Hog1SignalingData(IDS(i,conc)); data.tt(1)
        [time,Im,Is] = get_ImIs(data,-2,25+3,'Vol'); size(Im)
        shadedErrorBar(time,Im0-Im,Is,{'LineWidth', lw, 'color',cmap{ii+2}(count,:)}, 0.2); 
        box on; xlim([-2.5 32.5]); ylim([-.02 .22]); yticks([0:.1:1]); xticks([0:25:50]); xticklabels([]); yticklabels([]); 
        % if ii>=1; yticklabels([]); end
        count = count + 1; 
    end
    ii=ii+1;
end    

%% volume sum polynomials 
% t1, t2, t3, t5, t7
IDS = [[15:19]; [20:24]; [25:29]; [30:34]; [35:39]]; 

FM = [.1 .2 .4 .6 .8]; % final molarities
polynomials = [1 2 3 5 7]; % polynomial orders 

sumVOLm = NaN(numel(Im0), 5); 
sumVOLs = NaN(numel(Im0), 5); 

for conc=1:numel(FM) % concentrations
    a(conc) = (1+FM(conc)-0)^(1/25);
    count = 1; 
    for i = 1:numel(polynomials) % polynomials 
        data = Hog1SignalingData(IDS(i,conc)); data.tt(1)
        [time,Im,Is] = get_ImIs(data,-2,25+3,'Vol'); 
        VVm=Im0-Im; VVs=Is; 
        sumVOLm(:,conc) = nansum([sumVOLm(:,conc) ((25*log(a(conc))).^polynomials(i))*VVm./(FM(conc)*factorial(polynomials(i)))],2);        
        sumVOLs(:,conc) = nansum([sumVOLs(:,conc) ((25*log(a(conc))).^polynomials(i))*VVs./(FM(conc)*factorial(polynomials(i)))],2);  
    end    
end

hh=subplotHJ(3,3,3+3,dy,dx); cla(hh); hold on; %grid on
ii=1;
for i=4:5
	shadedErrorBar(time,sumVOLm(:,i),sumVOLs(:,i),{'LineWidth', lw, 'color',cmap{2}(ii,:)}, 0.2); 
    ii=ii+1; 
end
box on; xlim([-2.5 32.5]); ylim([-.02 .22]); yticks([0:.1:1]); xticks([0:25:50]); xticklabels([]); yticklabels([]); 
%if ii>=1; yticklabels([]); end

%% volume linear and exponential stimuli
IDS = [29 31; 95 97]; % t1 0.6M 0.8M; exp 0.6M 0.8M 

ii=1; 
for i=[1 2]
    hh=subplotHJ(3,3,6+ii,dy,dx); cla(hh); hold on; %grid on
    count = 1; 
    for exp = [1 2] %(2:5)
        data = DataSet1{IDS(i,exp)};
        shadedErrorBar(data.tt,data.mHog1nuc,data.sHog1nuc,{'LineWidth', lw, 'color',cmap{ii+2}(count,:)}, 0.2); 
        box on; xlim([-2.5 32.5]); ylim([-.02 1.22]); yticks([0:.5:1]); xticks([0:25:50]); xticklabels([]); yticklabels([]); 
%          if ii>=1; yticklabels([]); end
        count = count + 1; 
    end
    ii=ii+1;
end    

%% Hog1 sum polynomials 

IDs = [2 3:2:9]; 
polk = [3 4 5 6 7]; 
FM = [.1 .2 .4 .6 .8]; % final molarities

[time,Im0,Is0] = get_ImIs(DataSet1{5},-2,32,'Hog1'); 

sumHog1m = NaN(numel(Im0), 5); 
sumHog1s = NaN(numel(Im0), 5); 

for conc=1:5 
    a(conc) = (1+FM(conc)-0)^(1/25);
    count = 1; 
    for polynomials = [1 2 3 5 7] 
        exp = 11*(polk(count)-1)+IDs(conc);
        data = DataSet1{exp}; 
        [time,Im,Is] = get_ImIs(data,-2,32,'Hog1'); 
        sumHog1m(:,conc) = nansum([sumHog1m(:,conc) ((25*log(a(conc))).^polynomials)*Im./(FM(conc)*factorial(polynomials))],2);        
        sumHog1s(:,conc) = nansum([sumHog1s(:,conc) ((25*log(a(conc))).^polynomials)*Is./(FM(conc)*factorial(polynomials))],2);  
        count  = count + 1;
    end    
end

hh=subplotHJ(3,3,6+3,dy,dx); cla(hh); hold on; %grid on
ii=1;
for i=4:5
	shadedErrorBar(time,sumHog1m(:,i),sumHog1s(:,i),{'LineWidth', lw, 'color',cmap{2}(ii,:)}, 0.2); 
    ii=ii+1; 
end
box on; xlim([-2.5 32.5]); ylim([-.02 1.22]); xticks([0:25:50]); xticklabels([]); yticklabels([]); yticks([0:.5:1])
%% save figure 
set(findall(gcf,'-property','FontSize'),'FontSize',6, 'defaultTextFontSize',6, 'FontName', 'Helvetica')
figname = [dir_name, '/HOGt1expsum'];  
print(figname,'-depsc', '-r600');

%% plot pulsatile staircase, T = 16 min
path0 = '/Users/hosseinjashnsaz/Dropbox (VU Basic Sciences)/Hossein Jashnsaz/Experiments/Diverse_kinetics_Hog1YFP_TimeLapse2/#7 data structure/'; 
DATA = load([path0, 'Hog1SignalingData_pulse_staircase2.mat']);
Hog1SignalingData = DATA.Hog1SignalingData2;
% pulse full, pulse 2nd, pulse 3rd, stair full, stair 2nd, stair 3rd 

close; figure();clf; set(gcf, 'Units', 'centimeters', 'Position', [0 0 5.8 12], 'PaperUnits', 'centimeters', 'PaperSize', [5.8 12]); 
pulsestair_col = {[0 0 1], hex2rgb('090707')}; lww = {2,1}; lw=1; 
ctr_col = .6*[1 1 1]; %ctrl

% Pulsatile 
% 2min-16min-16min-16min-16min
% 0M - 0.2M - 0.4M - 0.6M-0.8 

% Staircase 
% 2min-8min8min-8min8min-8min8min
% 0M - 0.2M-0M - 0.2M-0M - 0.2M-0M
NaCls = {[0 0 0.2 0.2 0 0 0.2 0.2 0 0 .2 .2 0 0], [0 0 0.2 0.2 0.4 0.4 0.6 0.6 .8 .8 .8]}; 
times = {[-2 0 0 8 8 8*2 8*2 8*3 8*3 8*4 8*4 8*5 8*5 8*6], [-2 0 0 16 16 16*2 16*2 16*3 16*3 16*4 16*4]}; 

% stimuli 
hh=subplotHJ(4,1,1,dy,dx); cla(hh); hold on; %grid on
plot([-2 58],[0 0], 'color',ctr_col, 'LineWidth', 2*lw); 
for i=1:2
    plot(times{i},NaCls{i}, 'color', pulsestair_col{i}, 'LineWidth', lww{i}); 
end
box on; xlim([-2 16*3+2]); ylim([-0.01 .61]);
xticks([0:8:16*4]); yticks([0:.2:.8]); xticklabels([]); yticklabels([]); 

% volume change
hh=subplotHJ(4,1,4,dy,dx); cla(hh); hold on; %grid on

% control, 0M
[time0,Im0,Is0] = get_ImIs(data48,-1,49,'Vol'); 
shadedErrorBar(time0-2,Im0,Is0,{'LineWidth', lw, 'color',ctr_col}, 0.2); 

% pulsatile, 123
data=Hog1SignalingData(1); data.tt([1 end])
[time1,Im1,Is1] = get_ImIs(data,-3,47,'Vol');
shadedErrorBar(time1,Im1,Is1,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 
[MINM16,INDX16] = min(abs(time1-16)); [MINM32,INDX32] = min(abs(time1-2*16));
% pulsatile, pulse 2
data=Hog1SignalingData(2); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Vol');
shadedErrorBar(16+time,Im1(INDX16)+Im-Im(1),Is,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 
% pulsatile, pulse 3
data=Hog1SignalingData(3); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Vol');
shadedErrorBar(2*16+time,Im1(INDX32)+Im-Im(1),Is,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 

% staircase, 123
data=Hog1SignalingData(4); data.tt([1 end])
[time2,Im2,Is2] = get_ImIs(data,-3,47,'Vol'); 
shadedErrorBar(time2,Im2,Is2,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2);
[MINM16,INDX16] = min(abs(time2-16)); [MINM32,INDX32] = min(abs(time2-2*16));
% staircase, pulse 2
data=Hog1SignalingData(5); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Vol');
shadedErrorBar(16+time,Im2(INDX16)+Im-Im(1),Is,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2); 
% staircase, pulse 3
data=Hog1SignalingData(6); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Vol');
shadedErrorBar(2*16+time,Im2(INDX32)+Im-Im(1),Is,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2); 

box on; xlim([-2 16*3+2]); %ylim([-0.01 .61]);
xticks([0:8:16*4]); yticks([0:.2:.8]); xticklabels([]); yticklabels([]); 

% volume shrink
hh=subplotHJ(4,1,2,dy,dx); cla(hh); hold on; %grid on
shadedErrorBar(time0,Im0-Im0,Is0,{'LineWidth', lw, 'color',ctr_col}, 0.2); 
shadedErrorBar(time1,Im0-Im1,Is1,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 
shadedErrorBar(time2,Im0-Im2,Is2,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2); 
box on; xlim([-2 16*3+2]); ylim([-0.05 .21]);
xticks([0:8:16*4]); yticks([0:.1:.8]); xticklabels([]); yticklabels([]); 

% Hog1nuc
hh=subplotHJ(4,1,3,dy,dx); cla(hh); hold on; %grid on
[time,Im0,Is0] = get_ImIs(data48,-1,49,'Hog1'); 
shadedErrorBar(time-2,Im0,Is0,{'LineWidth', lw, 'color',ctr_col}, 0.2); 

% pulsatile, pulse 1
data=Hog1SignalingData(1); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Hog1'); mIm1=max(Im); 
shadedErrorBar(time,Im,Is,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 
% pulsatile, pulse 2
data=Hog1SignalingData(2); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Hog1');
shadedErrorBar(16+time,Im,Is,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 
% pulsatile, pulse 3
data=Hog1SignalingData(3); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Hog1');
shadedErrorBar(2*16+time,Im,Is,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 

% staircase, pulse 1
data=Hog1SignalingData(4); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Hog1'); mIm2=max(Im); mIm=max([mIm1 mIm2]); 
shadedErrorBar(time,Im,Is,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2); 
% staircase, pulse 2
data=Hog1SignalingData(5); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Hog1');
shadedErrorBar(16+time,Im,Is,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2); 
% staircase, pulse 3
data=Hog1SignalingData(6); data.tt([1 end])
[time,Im,Is] = get_ImIs(data,-3,14,'Hog1'); Is(Is>.5)=NaN; 
shadedErrorBar(2*16+time,Im,Is,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2); 

box on; xlim([-2 16*3+2]); ylim([-0.01 .71]);
xticks([0:8:16*4]); yticks([0:mIm/2:2]); xticklabels([]); yticklabels([]); 

set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica')
figname = [dir_name, '/pulse_stair'];  
print(figname,'-depsc', '-r600');

%% Figure S5  
close; figure();clf; set(gcf, 'Units', 'centimeters', 'Position', [0 0 20.5 14], 'PaperUnits', 'centimeters', 'PaperSize', [20.5 14]); 
dx = 0.03; dy = 0.05; lw=1; lww={3,2}; 

hh=subplotHJ(3,1,1,dy,dx); cla(hh); hold on; %grid on
for i=1:2
    plot(times{i},NaCls{i}, 'color', pulsestair_col{i}, 'LineWidth', lww{i}); 
end
box on; xlim([-1.1 16*3+1]); ylim([-0.01 .61]);
xticks([0:8:16*4]); yticks([0:.2:.8]); % xticklabels([]); yticklabels([]); 

% volume change
for i=1:3
hh=subplotHJ(3,3,3+i,dy,dx); cla(hh); hold on; %grid on
% control, 0M
% [time0,Im0,Is0] = get_ImIs(data48,-1+16*(i-1),16+16*(i-1),'Vol'); 
% shadedErrorBar(time0-2-16*(i-1),Im0,Is0,{'LineWidth', lw, 'color',ctr_col}, 0.2); 

% pulsatile, i
data=Hog1SignalingData(i); data.tt([1 end])
[time1,Im1,Is1] = get_ImIs(data,-3,14,'Vol');
shadedErrorBar(time1,Im1,Is1,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 

% staircase, i
data=Hog1SignalingData(3+i); data.tt([1 end])
[time2,Im2,Is2] = get_ImIs(data,-3,14,'Vol'); 
shadedErrorBar(time2,Im2,Is2,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2);

box on; xlim([-1.2 14.2-2]); ylim([0.89 1.11]);
xticks([0:8:16]); yticks([.9:.1:1.1]); xticklabels(16*(i-1)+[0:8:16]); % yticklabels([]); 
end

% Hog1nuc
for i=1:3
hh=subplotHJ(3,3,6+i,dy,dx); cla(hh); hold on; %grid on
% control, 0M
% [time0,Im0,Is0] = get_ImIs(data48,-1+16*(i-1),16+16*(i-1),'Hog1'); 
% shadedErrorBar(time0-2-16*(i-1),Im0,Is0,{'LineWidth', lw, 'color',ctr_col}, 0.2); 

% pulsatile, i
data=Hog1SignalingData(i); data.tt([1 end])
[time1,Im1,Is1] = get_ImIs(data,-3,14,'Hog1');
shadedErrorBar(time1,Im1,Is1,{'LineWidth', lw, 'color',pulsestair_col{1}}, 0.2); 

% staircase, i
data=Hog1SignalingData(3+i); data.tt([1 end])
[time2,Im2,Is2] = get_ImIs(data,-3,14,'Hog1'); 
shadedErrorBar(time2,Im2,Is2,{'LineWidth', lw, 'color',pulsestair_col{2}}, 0.2);

box on; xlim([-1.2 14.2-2]); ylim([-0.01 .71]); yticks([0:mIm/2:2]); 
xticks([0:8:16]); xticklabels(16*(i-1)+[0:8:16]); yticklabels([]); 
end

set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica')
figname = [dir_name, '/Figure_S5_pulse_stair'];  
print(figname,'-depsc', '-r600');

%% plot staircase 2
figure();clf; set(gcf, 'Units', 'centimeters', 'Position', [0 0 .3*13.6 6], 'PaperUnits', 'centimeters', 'PaperSize', [.3*13.6 6]); 
stair_col = hex2rgb('090707'); %step orange
ctr_col = .6*[1 1 1]; %ctrl

Experiments = {'0.30M25min 0.30M25min 0.90M25min', '0.20M25min 0.40M25min 0.60M25min', '0.20M25min 0.40M25min 0.60M25min 0.80M25min'}; 
tmax = [77, 80, 100]; lw = 1; dx=1*dx;

ExperimentsID = {[44], [45], [46]}; 
NaCls = {[0 0 0.3 0.3 0.6 0.6 0.9 0.9 .9], [0 0 0.2 0.2 0.4 0.4 0.6 0.6 .6], [0 0 0.2 0.2 0.4 0.4 0.6 0.6 0.8 0.8 .8]}; 
times = {[-2 0 0 25 25 50 50 75 77], [-2 0 0 25 25 50 50 75 tmax(1)], [-2 0 0 25 25 50 50 75.5 75.5 100 100]}; 
stairs = [.3 .2 .2]; 
pls=[1 4 7 2 5 8 3 6 9]; 

ii = 1; 
for i=size(ExperimentsID,2)
    subplotHJ(2,1,1,dy,dx); hold on; %grid on
    
    plot([-2 58],[0 0], 'color',ctr_col, 'LineWidth', lw); 
    plot(times{i},NaCls{i}, 'color', stair_col, 'LineWidth', lw); 
%     box on; xlim([-2 tmax(i)+1]); ylim([-0.01 .81]);
    box on; xlim([-2 75]); ylim([-0.01 .61]);
    xticks([0:25:tmax(i)]); yticks([0:stairs(i):.8]); xticklabels([]); 
    
    ii = ii + 1; 
    subplotHJ(2,1,2,dy,dx); hold on; %grid on
    count = 1; 
    for exp = ExperimentsID{i}
        
        % control 
        data = DataSet1{100}; 
        shadedErrorBar(data.tt,data.mHog1nuc,data.sHog1nuc,{'LineWidth', lw, 'color',ctr_col}, 0.2); 
                
        data = DataSet3{7}; 
        shadedErrorBar(data.tt(1:length(data.mHog1nuc)),data.mHog1nuc,data.sHog1nuc,{'LineWidth', lw, 'color',stair_col}, 0.2); 
%         box on; xlim([-2 tmax(i)+1]); ylim([-.02 1.22]); xticks([0:25:tmax(i)]);
        box on; xlim([-2 75]); ylim([-.02 1.02]); xticks([0:25:tmax(i)]); yticks([0:.5:1])
        count = count + 1; 
    end
    ii = ii + 1; 
end    

set(findall(gcf,'-property','FontSize'),'FontSize',6, 'defaultTextFontSize',6, 'FontName', 'Helvetica')
figname = [dir_name, '/HOGstair2'];  
print(figname,'-depsc', '-r600');

