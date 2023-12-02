%% Growth rate (OD600 measuremnts) analysis upon diverse kinetics
clc; close all; clear all; 
dir_name = 'plots'; mkdir(dir_name); 
addpath('../../../data/matlab-general-functions/')

day='day7'; logs = 0; meanBRs=0; smooth_mBRs=0; fullexp = 0; 
OD600ticks=[-3:1:0]; 
OD600ticksLabels=sprintf('2^{%d}\n',OD600ticks); 

% set # of time points for  exp growth curve fitting
optionN.drop="no";
optionY.drop="yes"; optionY.howmany = 4; 
options{1}.options{1}=optionN; %% pre-stimuli 
options{2}.options{1}=optionY; % options{2}.options{2}=optionN; % 1st measurments  
options{3}.options{1}=optionY; % options{3}.options{2}=optionN; % 2nd measurments  
options{4}.options{1}=optionN;  % 3rd measurments  

%% stimuli functions
%control stimuli u=u0;
constant_stimuli_input = @(t, conc_start) (t-t + conc_start);

% step stimuli
sympref('HeavisideAtOrigin',1);
step_stimuli_profile = @(t, conc_start, conc_change) conc_start + conc_change * heaviside(t - .001); 
% stimuli profile as polynomial functions
polynomial_stimuli_profile = @(t, conc_start, k, duration, conc_change) max(conc_start,min([conc_start + (conc_change.*t.^k/(duration)^k),conc_start+conc_change]));
% polynomial_stimuli_profile = @(t, conc_start, k, duration, conc_change) conc_start + (conc_change.*t.^k/(duration)^k);
% stimuli profile as exponential functions
exponential_stimuli_profile = @(t, conc_start, duration, conc_change) max(conc_start, min([(conc_start -1 + ((1+conc_change)^(1/(duration)))^t),conc_start+conc_change]));

cmap = HJ_costum_colors(3,.35);  
nfig =0; DATA = []; nfig = 1; 
%% set directories
% % global data path
% % path0 = '/Users/hosseinjashnsaz/NeuertLab/Dropbox (VU Basic Sciences)/Hossein Jashnsaz/Experiments/'; 
% path0 = '/Users/hosseinjashnsaz/Dropbox (VU Basic Sciences)/Hossein Jashnsaz/Experiments/'; 
% subpath = 'Diverse_kinetics_Hog1YFP_TimeLapse/'; 
% Folder = '#9 Growth rates/';
% 
% path_to_data  = [path0, subpath, Folder];
% addpath(path_to_data)

% local data path
addpath('../../../data/cell-growth-data-OD600-measurements/')
%% load OD600 data day7 20200125 ctrl, 120min pulse, 60min pulse, 90min pulse, 30min pulse
OD600s=importdata(['HJ_20200125_WT_BY4741_data.csv']);             % load data 
OD600s = OD600s.data; 
t0(1) = 98; % min
EXP_ID{1}=[1:5];  

DATA = [DATA OD600s]; clear OD600s; 
%% load OD600 data day8 20200126 ctrl, 30min pulse 30min shift, 30min pulse 60min shift, 30min pulse 90min shift, t1-90min-3M-30min
OD600s=importdata(['HJ_20200126_WT_BY4741_data.csv']);             % load data 
OD600s = OD600s.data;
t0(2) = 153; 
EXP_ID{2}=[6:10]; 

DATA(:,size(DATA,2)+1) = NaN; 
if size(OD600s,1)>size(DATA,1)
    DATA(size(DATA,1)+1:size(OD600s,1),:) = NaN; 
else
    OD600s(size(OD600s,1)+1:size(DATA,1),:) = NaN; 
end

DATA = [DATA OD600s]; clear OD600s; 

%% load OD600 data day5 20191214 ctrl, exp CSM, pulse20, exp_20m, exp_30m
OD600s=importdata(['HJ_20191214_WT_BY4741_data.csv']);             % load data 
OD600s = OD600s.data; 
t0(3) = 105; 
EXP_ID{3}=[11:12];  

DATA(:,size(DATA,2)+1) = NaN; 
if size(OD600s,1)>size(DATA,1)
    DATA(size(DATA,1)+1:size(OD600s,1),:) = NaN; 
else
    OD600s(size(OD600s,1)+1:size(DATA,1),:) = NaN; 
end
DATA = [DATA OD600s(:,[1:5 21:24])]; clear OD600s; %get cntrl and exp-30min 

OD600s  = DATA; 
stimuli_start_time = nanmax(t0); % min

if fullexp; t_max = nanmax(nanmax(OD600s)); else; t_max= stimuli_start_time+14*60;end; dt = .1;
time_points = [0:dt:t_max]; 

%% stimuli inputs
stimuli = zeros(20,size(time_points,2)); 

PULSE1_120min3M = find((stimuli_start_time)<=time_points & time_points<=(stimuli_start_time+120) | (stimuli_start_time+120+300)<=time_points & time_points<=(stimuli_start_time+120+300+120)); 
PULSE2_090min3M = find((stimuli_start_time)<=time_points & time_points<=(stimuli_start_time+90) | (stimuli_start_time+120+300)<=time_points & time_points<=(stimuli_start_time+120+300+90)); 
PULSE3_060min3M = find((stimuli_start_time)<=time_points & time_points<=(stimuli_start_time+60) | (stimuli_start_time+120+300)<=time_points & time_points<=(stimuli_start_time+120+300+60)); 
PULSE4_030min3M = find((stimuli_start_time)<=time_points & time_points<=(stimuli_start_time+30) | (stimuli_start_time+120+300)<=time_points & time_points<=(stimuli_start_time+120+300+30)); 
IN_ID=1; 

stimuli(IN_ID,:)=0; IN_ID=IN_ID+1; 
stimuli(IN_ID,PULSE1_120min3M)=3; IN_ID=IN_ID+1; 
stimuli(IN_ID,PULSE2_090min3M)=3; IN_ID=IN_ID+1;
stimuli(IN_ID,PULSE3_060min3M)=3; IN_ID=IN_ID+1;
stimuli(IN_ID,PULSE4_030min3M)=3; IN_ID=IN_ID+1;


PULSE5_30min3M_30min = find((stimuli_start_time+30)<=time_points & time_points<=(stimuli_start_time+30+30) | (stimuli_start_time+120+300+30)<=time_points & time_points<=(stimuli_start_time+120+300+30+30)); 
PULSE6_30min3M_60min = find((stimuli_start_time+60)<=time_points & time_points<=(stimuli_start_time+60+30) | (stimuli_start_time+120+300+60)<=time_points & time_points<=(stimuli_start_time+120+300+60+30)); 
PULSE7_30min3M_90min = find((stimuli_start_time+90)<=time_points & time_points<=(stimuli_start_time+90+30) | (stimuli_start_time+120+300+90)<=time_points & time_points<=(stimuli_start_time+120+300+90+30)); 
t1_90min3M = find(stimuli_start_time<=time_points & time_points<=(stimuli_start_time+90) | (stimuli_start_time+120+300)<=time_points & time_points<=(stimuli_start_time+120+300+90)); %  2, 310
t1EXT30min3M = find((stimuli_start_time+90)<=time_points & time_points<=(stimuli_start_time+90+30) | (stimuli_start_time+120+300+90)<=time_points & time_points<=(stimuli_start_time+120+300+90+30)); 

stimuli_input = @(t) polynomial_stimuli_profile(t, 0, 1, 90, 3); 
stimuli_TPs = 0:dt:90;

for t=1:length(stimuli_TPs)
    t1_stimuli(t) = stimuli_input(stimuli_TPs(t)); 
end

stimuli(IN_ID,:)=0; IN_ID=IN_ID+1; 
stimuli(IN_ID,PULSE5_30min3M_30min)=3; IN_ID=IN_ID+1; 
stimuli(IN_ID,PULSE6_30min3M_60min)=3; IN_ID=IN_ID+1; 
stimuli(IN_ID,PULSE7_30min3M_90min)=3; IN_ID=IN_ID+1; 
stimuli(IN_ID,t1_90min3M)=[t1_stimuli  t1_stimuli]; stimuli(IN_ID,t1EXT30min3M)=3; IN_ID=IN_ID+1; 

% te_90min3M = find(stimuli_start_time<=time_points & time_points<=(stimuli_start_time+90) | (stimuli_start_time+90+20+320)<=time_points & time_points<=(stimuli_start_time+90+20+320+90)); 
te_90min3M = find(stimuli_start_time<=time_points & time_points<=(stimuli_start_time+90) | (stimuli_start_time+120+300)<=time_points & time_points<=(stimuli_start_time+120+300+90)); 
teEXT30min = find(stimuli_start_time+90<time_points & time_points<=(stimuli_start_time+90+30) | (stimuli_start_time+120+300+90)<time_points & time_points<=(stimuli_start_time+120+300+90+30)); 
stimuli_input = @(t) exponential_stimuli_profile(t, 0, 90, 3); 
for t=1:length(stimuli_TPs)
    te_stimuli(t) = stimuli_input(stimuli_TPs(t)); 
end

stimuli(IN_ID,:)=0; IN_ID=IN_ID+1; 
stimuli(IN_ID,te_90min3M)=[te_stimuli  te_stimuli]; stimuli(IN_ID,teEXT30min)=3; IN_ID+1; 
%% figure
nfig=nfig+1;figure(nfig); set(gcf,'defaultLineLineWidth',1); set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 1]);
set(gcf, 'Units', 'centimeters', 'Position', [0 0 21 22], 'PaperUnits', 'centimeters', 'PaperSize', [21 22]); 

dx = .05; dy = .01;  
% cols=[7,6,4,9,8, 7,6,4,9,8, 7,6,4,9,8, 7,6,4,9,8];
cols = [12 13 13 12 13 13 12 11 9]; fits_col=16;  
ii = 1; 
for i=[1 6 11]
    subplotHJ(10,3,ii,dy,dx); grid on; hold on; box on;       
%     yyaxis left % Growth OD600
    times = (i-1)*5+1; 
    for eee=1:3 
        if ismember(i,EXP_ID{eee})==1;  shift_times = stimuli_start_time-t0(eee);end
    end
    ODs = []; 
    for BRs=1:3
        BR_Measurements_shifts = (BRs-1)+.5; TP_min = shift_times + OD600s(:,times)+BR_Measurements_shifts; TP_h = TP_min/60; 
        ODs = [ODs OD600s(:,times+BRs)]; 
        ss=scatter(TP_h, log2(OD600s(:,times+BRs)),60); ss.MarkerFaceColor = cmap{1}(BRs,:); ss.MarkerEdgeColor = 'none'; ss.MarkerFaceAlpha = .7; 
        [timepoints, OD600] = section_each_treatment(TP_min',OD600s(:,times+BRs)'); 
        for k=1%:size(OD600,2)
             tt_k = [timepoints{k}]'; OD_k = [OD600{k}]'; %OD_k = smooth(tt_k,OD_k,0.3,'rloess');
             for n_fits=1:100
                [t00,f,CI_f,mean_T,neg_CI,pos_CI] = get_fits(tt_k, OD_k, optionN);
                TT(n_fits) = mean_T; 
                if (k==1); start_ODs(BRs,n_fits) = f(stimuli_start_time-tt_k(1)); end
             end
             tps = [linspace(t00,tt_k(end),100)]'; 
             plot(tps/60,log2(f(tps-t00)), '-', 'color', [cmap{fits_col}(BRs,:) .6]);
             GROWTHR0M(ii,k,BRs)=nanmean(TT); 
        end
        
    end
    xlim([-.5 4.5]); xticks([0:.5:4]); xticklabels([]); 
    if ii==1;ylabel('OD600');end; if logs; ylim([.08 1.1]); set(gca, 'Yscale', 'log'); else; ylim(log2([.18 1.1])); yticks([OD600ticks]); yticklabels({OD600ticksLabels}); end
    ii=ii+1; 
end

% selects = [2 5 9 10 12]; % 2h .5h .5 t1 exp
selects = [2:5 7:10 12];
sections_pos = [stimuli_start_time/60-1.5, stimuli_start_time/60+4.5, stimuli_start_time/60+4.5+7]; 
OD600ticks=[-3:1:-1]; 
OD600ticksLabels=sprintf('2^{%d}\n',OD600ticks); 
ii = 1; 
for i=selects
    set(gcf,'defaultAxesColorOrder',[0 0 0; cmap{cols(ii)}(1,:)]);   
    subplotHJ(10,1,1+ii,dy,dx); grid on; hold on; 
    yyaxis left % Growth OD600
    times = (i-1)*5+1; 
    for eee=1:3 
        if ismember(i,EXP_ID{eee})==1; eeee=eee; shift_times = stimuli_start_time-t0(eee);end
    end
    ODs = []; 
    for BRs=1:3
        BR_Measurements_shifts = (BRs-1)+.5; TP_min = shift_times + OD600s(:,times)+BR_Measurements_shifts; TP_h = TP_min/60; 
        ODs = [ODs OD600s(:,times+BRs)]; 
        ss=scatter(TP_h, log2(OD600s(:,times+BRs)),60); ss.MarkerFaceColor = cmap{1}(BRs,:); ss.MarkerEdgeColor = 'none'; ss.MarkerFaceAlpha = .7; 
        [timepoints, OD600] = section_each_treatment(TP_min',OD600s(:,times+BRs)'); 
        for k=1:size(OD600,2)
             tt_k = [timepoints{k}]'; OD_k = [OD600{k}]'; %OD_k = smooth(tt_k,OD_k,0.3,'rloess');
             for kk=1:length(options{k}.options)
                 option = options{k}.options{kk}; 
%                  if (k==1 & ismember(ii,[5:8])==1); option=optionY; option.howmany=3; end
                 for n_fits=1:100
                    [t00,f,CI_f,mean_T,neg_CI,pos_CI] = get_fits(tt_k, OD_k, option);
                    TT(n_fits) = mean_T; 
                    if (k==1); start_ODs(BRs,n_fits) = f(stimuli_start_time-t00); end
                 end
                 tps = [linspace(t00,tt_k(end),100)]'; 
                 plot(tps/60,log2(f(tps-t00)), '-', 'color', [cmap{fits_col+kk-1}(BRs,:) .6]);
                 if nanmean(TT)<0; TT=NaN*TT; end; 
                 GrowthRates{ii,BRs}.sections(k)=nanmean(TT); 
                 GrowthRatesSTD{ii,BRs}.sections(k)=nanstd(TT);
                 GROWTHR(ii,k,BRs,kk)=nanmean(TT); 
             end
        end     
    end
        
    % means of BRs
    TTT=[]; 
    if meanBRs
    ss=scatter(TP_h, nanmean(ODs,2),60,'s'); ss.MarkerFaceColor = cmap{15}(1,:); ss.MarkerEdgeColor = 'none'; ss.MarkerFaceAlpha = .7; 
    [timepoints, OD600] = section_each_treatment(TP_min',nanmean(ODs,2)'); 
        for k=1:size(OD600,2)
             if k==1; option=optionN; else; option=optionY; end
             tt_k = [timepoints{k}]'; OD_k = [OD600{k}]'; 
             for n_fits=1:10
                [t00,f,CI_f,mean_T,neg_CI,pos_CI] = get_fits(tt_k, OD_k, option);
                tps = [linspace(t00,tt_k(end),100)]'; plot(tps/60,f(tps-t00), '-', 'color', [cmap{fits_col}(1,:) .6]);
                TT(n_fits) = mean_T; 
                
                if smooth_mBRs
                    smooth_OD_k = smooth(OD_k);
                    s=scatter(tt_k/60,smooth_OD_k,50,'s'); s.MarkerFaceColor = 'none'; s.MarkerEdgeColor = cmap{fits_col}(1,:); s.MarkerFaceAlpha = .3; 
                    [t00, f,~,mean_T,~,~] = get_fits(tt_k, smooth_OD_k, options);
                    plot(tps/60,f(tps-t00), '-', 'color', [cmap{16}(1,:) .3]);
                    TTT(n_fits) = mean_T; 
                end

             end
             GrowthRates{ii,BRs+1}.sections(k)=nanmean(TT); GrowthRates{ii,BRs+2}.sections(k)=nanmean(TTT);
             GrowthRatesSTD{ii,BRs+1}.sections(k)=nanstd(TT); GrowthRatesSTD{ii,BRs+2}.sections(k)=nanstd(TTT); 
        end        
    end
    
    textY=log2([.55 .4]); if ii==7; textY=log2([.20 .20]);end 
    for k=1:3
        if k==1;           
            grt=text(sections_pos(k), log2(.55), [num2str(round(nanmean(GROWTHR0M(eeee,k,:)),0)),'+/-',  num2str(max([round(nanstd(GROWTHR0M(eeee,k,:)),0) 1]))], 'color', cmap{fits_col}(1,:)); grt.HorizontalAlignment='c';
        else
            grt=text(sections_pos(k), textY(1), [num2str(round(nanmean(GROWTHR(ii,k,:,1)),0)),'+/-',  num2str(round(nanstd(GROWTHR(ii,k,:,1)),0))],'color', cmap{fits_col}(1,:)); grt.HorizontalAlignment='c'; 
%             grt=text(sections_pos(k), textY(2), [num2str(round(nanmean(GROWTHR(ii,k,:,2)),0)),'+/-',  num2str(round(nanstd(GROWTHR(ii,k,:,2)),0))],'color', cmap{fits_col+1}(1,:)); grt.HorizontalAlignment='c'; 
        end
    end
%     shadedErrorBar(OD600s(:,times)/60, nanmean(ODs,2),nanstd(ODs,0,2), {'-','LineWidth', 2, 'color', cmap{1}(1,:)},.15);
    xlim([-.5 ceil(t_max/60)]); xticks([(stimuli_start_time/60)-4:.5:20]); xticklabels([]); 
    ylabel('OD600'); if logs; set(gca, 'Yscale', 'log'); else; ylim(log2([.07 .8])); yticklabels({OD600ticksLabels}); end
    save('GROWTHR0M','GROWTHR0M');
    yyaxis right; % Stimuli INPUT (M)
    ylim([0, 3.2]); yticks([0 3]); ylabel('NaCl (M)'); %,'color','b'
    plot(time_points/60, stimuli(i,:), '-','color', cmap{cols(ii)}(1,:)); % cmap{cols(i)}(1,:)
    h=fill(time_points/60, stimuli(i,:), cmap{cols(ii)}(1,:));  h.EdgeColor = 'non'; h.FaceAlpha=.15; 
    set(gca,'GridLineStyle',':'); 
    box on; 
    
    yyaxis left; % Stimuli INPUT (M)
    start_OD = round(nanmean(nanmean(start_ODs,1)),3); txt=text(stimuli_start_time/60,log2(.15),num2str(start_OD)); txt.HorizontalAlignment='c'; 
    ii = ii +1; 
end
% xlabel('time (h)');
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

figname = [dir_name, '/Growth2','_LOG', num2str(logs), '_mu',num2str(meanBRs), '_S',num2str(smooth_mBRs), '_F', num2str(fullexp)];  
% print(figname,'-depsc', '-r600'); 
print(figname,'-dpng'); 

%% doubling times
nfig=nfig+1;figure(nfig); set(gcf,'defaultLineLineWidth',1);
set(gcf, 'Units', 'centimeters', 'Position', [0 0 4 21], 'PaperUnits', 'centimeters', 'PaperSize', [4 21]);
dx = .2; dy = .02;  

DoublingTimesExps=NaN(10,3,10); 
for i=1:ii-1
    DoublingTimes=NaN(3,10); % (BR, sec)
    subplotHJ(ii-1,1,i,dy,dx); hold on; grid on; 
    for BRs=1:3
        nn=length(GrowthRates{i,BRs}.sections);
        errorbar([1:length(GrowthRates{i,BRs}.sections)]+.2*(BRs-2),GrowthRates{i,BRs}.sections,GrowthRatesSTD{i,BRs}.sections,'.','color', cmap{1}(BRs,:)); 
        ss=scatter([1:length(GrowthRates{i,BRs}.sections)]+.2*(BRs-2), GrowthRates{i,BRs}.sections); 
        ss.MarkerFaceColor = cmap{1}(BRs,:); ss.MarkerEdgeColor = cmap{1}(BRs,:); 
        DoublingTimes(BRs,[1:nn]) = GrowthRates{i,BRs}.sections; %([1:min(3,length(GrowthRates{i,BRs}.sections))])
    end
    
    % DT estimated from means of BRs
    if meanBRs
       ss=scatter([1:length(GrowthRates{i,BRs+1}.sections)]+.4, GrowthRates{i,BRs+1}.sections,50, 'd'); 
       ss.MarkerEdgeColor = cmap{15}(BRs,:); ss.MarkerFaceColor = cmap{15}(BRs,:);  
       errorbar([1:length(GrowthRates{i,BRs+1}.sections)]+.4,GrowthRates{i,BRs+1}.sections,GrowthRatesSTD{i,BRs+1}.sections,'.','color', cmap{15}(BRs,:));
       
       s=scatter([1:length(GrowthRates{i,BRs+2}.sections)]+.5, GrowthRates{i,BRs+2}.sections,60, 'p'); 
       s.MarkerEdgeColor = cmap{16}(1,:); s.MarkerFaceColor = cmap{16}(1,:);  
       errorbar([1:length(GrowthRates{i,BRs+2}.sections)]+.5,GrowthRates{i,BRs+2}.sections,GrowthRatesSTD{i,BRs+2}.sections,'.','color', cmap{16}(1,:));
    end
    
    if fullexp; nn=length(GrowthRates{i,1}.sections); else; nn=3; end
    boxplot(DoublingTimes); 
    xlim([.4 nn+.6]); xticks([1:nn]); ylim([90 310]); yticks([100 200 300]);  set(gca,'YScale', 'log'); %ylabel('Doubling Time (min)'); %ylim([100 3600]); yticks([100 500 1000 3000]);grid on; set(gca,'YScale','log')
    DoublingTimesExps(i,:,:)=DoublingTimes; 
end
% xlabel('section');
save('DoublingTimesExps','DoublingTimesExps'); 

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'defaultTextFontSize',12, 'FontName', 'Helvetica')

figname = [dir_name, '/DoublingTime2','_LOG', num2str(logs), '_mu',num2str(meanBRs), '_S',num2str(smooth_mBRs), '_F', num2str(fullexp)];  
print(figname,'-dpng'); 

%% doubling rate
nfig=nfig+1;figure(nfig); set(gcf,'defaultLineLineWidth',1);
set(gcf, 'Units', 'centimeters', 'Position', [0 0 4 21], 'PaperUnits', 'centimeters', 'PaperSize', [4 21]);
dx = .2; dy = .02;  
PRESTIMULI = GROWTHR0M(:); 

% control, steps, linear, exp
mPRESTIMULI = nanmean(PRESTIMULI); 
controls = PRESTIMULI/mPRESTIMULI;
controls=1./controls; 

nDoublingTimesExps = NaN(9,3,2);  % normalize doubling times to control
selects = 1:9; % [1 4 7 8 9]; 
for i=1:length(selects)
    for j=1:3
        nDoublingTimesExps(i,:,j) = DoublingTimesExps(selects(i),:,j)/mPRESTIMULI; 
    end
end
DoublingRates =1./nDoublingTimesExps; 

for i=1:ii-1
    subplotHJ(ii-1,1,i,dy,dx); hold on; grid on; 
    for BRs=1:3
%         errorbar([1:length(DoublingRates(i,BRs,:))]+.2*(BRs-2),DoublingRates(i,BRs,:),GrowthRatesSTD{i,BRs}.sections,'.','color', cmap{1}(BRs,:)); 
        ss=scatter([1:length(DoublingRates(i,BRs,:))]+.2*(BRs-2), reshape(DoublingRates(i,BRs,:),[1 3])); 
        ss.MarkerFaceColor = cmap{1}(BRs,:); ss.MarkerEdgeColor = cmap{1}(BRs,:); 
    end
    thisDoublingRates=reshape(DoublingRates(i,:,:),[3 3]); % (BR, section)
    boxplot(thisDoublingRates); 
    xlim([.5 3.5]); xticks([1.5:1:3.5]); xticklabels([]); ylim([-.05 1.05]); yticks([0:.2:1]);  set(gca,'YScale', 'lin'); %ylabel('Doubling Time (min)'); %ylim([100 3600]); yticks([100 500 1000 3000]);grid on; set(gca,'YScale','log')
end

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'defaultTextFontSize',12, 'FontName', 'Helvetica')

figname = [dir_name, '/DoublingRate2','_LOG', num2str(logs), '_mu',num2str(meanBRs), '_S',num2str(smooth_mBRs), '_F', num2str(fullexp)];  
print(figname,'-dpng'); 

% %% doubling rate
% nfig=nfig+1;figure(nfig); set(gcf,'defaultLineLineWidth',1);
% set(gcf, 'Units', 'centimeters', 'Position', [0 0 4 21], 'PaperUnits', 'centimeters', 'PaperSize', [4 21]);
% dx = .2; dy = .02;  
% PRESTIMULI = GROWTHR0M(:); 
% 
% % control, steps, linear, exp
% % mPRESTIMULI = nanmean(PRESTIMULI); 
% controls = 1./PRESTIMULI;
% mPRESTIMULI = nanmean(controls);
% controls=controls/mPRESTIMULI; 
% 
% nDoublingTimesExps = NaN(9,3,2);  % normalize doubling times to control
% selects = 1:9; % [1 4 7 8 9]; 
% for i=1:length(selects)
%     for j=1:3
%         nDoublingTimesExps(i,:,j) = 1./DoublingTimesExps(selects(i),:,j); 
%     end
% end
% DoublingRates =nDoublingTimesExps/mPRESTIMULI; 
% 
% for i=1:ii-1
%     subplotHJ(ii-1,1,i,dy,dx); hold on; grid on; 
%     for BRs=1:3
% %         errorbar([1:length(DoublingRates(i,BRs,:))]+.2*(BRs-2),DoublingRates(i,BRs,:),GrowthRatesSTD{i,BRs}.sections,'.','color', cmap{1}(BRs,:)); 
%         ss=scatter([1:length(DoublingRates(i,BRs,:))]+.2*(BRs-2), reshape(DoublingRates(i,BRs,:),[1 3])); 
%         ss.MarkerFaceColor = cmap{1}(BRs,:); ss.MarkerEdgeColor = cmap{1}(BRs,:); 
%     end
%     thisDoublingRates=reshape(DoublingRates(i,:,:),[3 3]); % (BR, section)
%     boxplot(thisDoublingRates); 
%     xlim([.5 3.5]); xticks([1.5:1:3.5]); xticklabels([]); ylim([-.05 1.05]); yticks([0:.2:1]);  set(gca,'YScale', 'lin'); %ylabel('Doubling Time (min)'); %ylim([100 3600]); yticks([100 500 1000 3000]);grid on; set(gca,'YScale','log')
% end
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12, 'defaultTextFontSize',12, 'FontName', 'Helvetica')
% 
% figname = [dir_name, '/DoublingRate3','_LOG', num2str(logs), '_mu',num2str(meanBRs), '_S',num2str(smooth_mBRs), '_F', num2str(fullexp)];  
% print(figname,'-dpng'); 

%% plot all together
for j=1:3

nfig=nfig+1;figure(nfig); clf; set(gcf, 'Units', 'centimeters', 'Position', [0 0 6 3], 'PaperUnits', 'centimeters', 'PaperSize', [6 3]); 
dx = 0.005; dy = 0.005; ssz=20; 
hold on

i=1; 
this_data = controls; 
bar(i,nanmean(this_data), 'FaceColor',cmap{4}(1,:),'EdgeColor','none','LineWidth',1.5);
errorbar(i, nanmean(this_data), nanstd(this_data),'LineWidth',2,'color', cmap{3}(1,:));

shifts = .15*[-1 0 1]; 
s=scatter(i+.2*randn(size(this_data)),this_data,ssz, 'o','LineWidth',.3); 
s.MarkerEdgeColor = cmap{4}(1,:); s.MarkerEdgeAlpha=.7; s.MarkerFaceColor = [cmap{4}(2,:)]; s.MarkerFaceAlpha=.7;  

selects=[NaN 1 2 3 4 5 6 7 8 9];
% selects=[NaN 1 2 3 4 5]; 

for i=2:length(selects)
    this_data = reshape(DoublingRates(selects(i),:,j),[1 3]);
    bar(i,nanmean(this_data), 'FaceColor',cmap{cols(i-1)}(1,:),'EdgeColor','none','LineWidth',1.5);
    errorbar(i, nanmean(this_data), nanstd(this_data),'LineWidth',2,'color', cmap{3}(1,:));
%     X=mvnrnd(mean(this_data),.5*std(this_data),100);
%     s=scatter(i+.1*randn(size(X)),X,2, 'o'); 
%     s.MarkerEdgeColor = cmap{i}(3,:); s.MarkerFaceColor = cmap{i}(3,:);%cmap{i}(3,:);  
    s=scatter(i+shifts,this_data,ssz, 'o'); 
    s.MarkerEdgeColor = cmap{4}(1,:); s.MarkerEdgeAlpha=.7; s.MarkerFaceColor = [cmap{4}(2,:)]; s.MarkerFaceAlpha=.7;  
end
hold on; box on; grid on; 
xlim([.5 length(selects)+.5]); xticks([1:length(selects)]); xticklabels([]); ylim([0 1.1]); yticks([0:.2:1]);
set(gca,'GridLineStyle',':'); 
% X1 = CONDs(4,:,1); X1=mvnrnd(mean(X1),std(X1),100);
% X2 = CONDs(5,:,1); X2=mvnrnd(mean(X2),std(X2),100);
% % [h,p] = kstest2(X1,X2)
% [h,p] = ttest2(X1,X2)

% % % X1 = 1+randn(1000,1); 
% % % X2 = 1.15+randn(1000,1); 
% % % [h,p] = kstest2(X1,X2)
% % % figure; boxplot([X1 X2])

set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

figname = [dir_name,'/growth_3_tp_', num2str(j)];  
print(figname,'-depsc', '-r600');

end

save('DATA_FINAL_3tp','controls','DoublingRates'); 
