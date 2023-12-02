warning('all','off')
% load Figure03
dir_name = 'plots'; mkdir(dir_name);
close all
addpath('../../../data/matlab-general-functions/')

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
addpath('../../../data/cell-volume-signaling-data/')
load Hog1SignalingData2
load HogDataSet 

%% define data
Hog1SignalingData=Hog1SignalingData2; 
DataSet1 = HogDataSet.DataSet1;

%% plot stimuli all polynomials
% close all
% figure(); set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 7.5], 'PaperUnits', 'centimeters', 'PaperSize', [20 7.5]); 
% dx = 0.035; dy = 0.06; lw=2; 
% 
% cmap = jet(8); 
% cmap = cmap([1:4 7:8],:); 
% cmap=[0 0 0; cmap]; 
% 
% for i=1:size(ExperimentsID,2)
%     subplotHJ(2,3,i,dy,dx); hold on; %grid on
% %     set(groot,'defaultAxesTickLabelInterpreter','latex');
%     text(8,.7, Experiments{i}, 'Interpreter','latex', 'HorizontalAlignment', 'c'); 
%     
%     count = 1;     
%     for exp = ExperimentsID{i}
%         plot(data{exp}.tt,data{exp}.stimulus_profile,'LineWidth', lw, 'color',cmap(count,:)); 
%         box on; xlim([-.5 50.5]); ylim([-.01 .81]); xticks([0:10:50]); yticks([0:.2:.8]); 
%         count = count + 1; 
%     end
% end
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')
% figname = [dir_name, '/NaCl'];  
% print(figname,'-depsc', '-r600');
%% figure handle
close; figure(1); set(gcf, 'Units', 'centimeters', 'Position', [0 0 21.5 6], 'PaperUnits', 'centimeters', 'PaperSize', [21.5 8]); 
dx = 0.05; dy = 0.08; cmap = jet(7); lw=2; cmap([3 5],:)=.95*cmap([3 5],:);

%% stimuli 
% data = stimuli_functions(); 
subplotHJ(1,4,1,dy,dx); hold on; grid on;
set(groot,'defaultAxesTickLabelInterpreter','tex');
count = 1; 
IDS  = [7:11:73];  
data = DataSet1{100}; 
plot(data.tt,-.1+data.stimulus_profile,'LineWidth', lw, 'color','k');
for i=IDS 
    data = DataSet1{i}; 
    timepoints = [data.tt(1):6/60:data.tt(end)]';% every 6 seconds (based on 25 minutes over 251 syringe pump phases)
    plot(timepoints,-.1+data.stimulus(timepoints),'LineWidth', lw, 'color',cmap(count,:)); 
    box on; xlim([-.5 50.5]); ylim([-.01 .81]); xticks([0:10:50]); yticks([0:.2:.8]); 
    count = count + 1; 
end
box on; xlim([-3.5 50.5]); ylim([-.01 .61]); xticks([0:10:50]); yticks([0:.1:.6]); 
set(gca,'GridLineStyle',':'); 
%% volume
Experiments = {'steps $ (t^{0}) $', 'root $ (\sqrt{t}) $', 'linears $ (t^{1}) $', 'quadratics $ (t^{2}) $', 'cubics $ (t^{3}) $', 'quintics $ (t^{5}) $', 'heptics $ (t^{7}) $', 'powers $ (a^t) $', 'sum polys'};
finalConcs = {[0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 1.0], [0.10 0.20:0.2:0.8], [[0.10 0.20:0.2:0.8]], ...
    [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]]}; 

concentrations = {'0.1M', '0.2M', '0.4M', '0.6M', '0.8M'}; 
ExperimentsID = {[1:3 6 8], [10:14], [15:19], [20:24], [25:29], [30:34], [35:39], [40:44]}; 

dd = 1;

for i=4%size(ExperimentsID,2)
    subplotHJ(1,4,3,dy,dx); hold on; grid on
    set(groot,'defaultAxesTickLabelInterpreter','tex');
%     text(25,1.25, concentrations{i}, 'Interpreter','latex', 'HorizontalAlignment', 'c'); 
    Hog1Volume.Vol = NaN(1000,8); 
    Hog1Volume.time = NaN(1000,8); 

    % control 
    exp = 48; 
%     n_BRs = size(Hog1SignalingData(exp).Sel,2); 
    [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
    shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color','k'}, 0.2); 
    Hog1Volume.Vol(1:numel(dp),1) = Hog1SignalingData(exp).Volm(dp);
    Hog1Volume.time(1:numel(dp),1) = Hog1SignalingData(exp).tt(dp); 
    Hog1Volume.TimeDelay(1) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 

    count = 1;     
    for j=[1:4 5 6:7]
        exp = ExperimentsID{j}(i); 
%         n_BRs = size(Hog1SignalingData(exp).Sel,2); 
        [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
        shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
        Hog1Volume.Vol(1:numel(dp),count+1) = Hog1SignalingData(exp).Volm(dp);
        Hog1Volume.time(1:numel(dp),count+1) = Hog1SignalingData(exp).tt(dp); 
        Hog1Volume.TimeDelay(count+1) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 
        count = count + 1; 
    end
    box on; xlim([-3.5 50.5]); ylim([.72 1.27]); xticks([0:10:50]);
    set(gca,'GridLineStyle',':'); 
end

%% Hog1 

    subplotHJ(1,4,4,dy,dx); hold on; grid on
    set(groot,'defaultAxesTickLabelInterpreter','tex');
%     text(25,.45, concentrations{i}, 'Interpreter','latex', 'HorizontalAlignment', 'c'); 
    Hog1Volume.Hog1 = NaN(1000,7);

    % control 
    data = DataSet1{100}; 
    shadedErrorBar(data.tt,data.mHog1nuc,data.sHog1nuc,{'LineWidth', lw, 'color','k'}, 0.2); 
    Hog1Volume.Hog1(1:numel(data.mHog1nuc),1) = data.mHog1nuc;
%     txt=text(48, 0.15, [num2str(0),'M']); txt.Color = 'k'; txt.FontWeight = 'b';
    count = 1; 
    for j=IDS
        data = DataSet1{j}; 
        shadedErrorBar(data.tt(1:length(data.mHog1nuc)),data.mHog1nuc,data.sHog1nuc,{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
        Hog1Volume.Hog1(1:numel(data.mHog1nuc),count+1) = data.mHog1nuc;
        count = count + 1; 
    end
    box on; xlim([-3.5 50.5]); ylim([-.02 1.48]); xticks([0:10:50]); 
    set(gca,'GridLineStyle',':'); 
%% save figs
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

figname = [dir_name, '/stimuli_cellvolume_Hog1'];  
print(figname,'-depsc', '-r600');
%% save data for this figure
% save('Figure03','-v7.3');


