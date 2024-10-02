warning('all','off')
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

%% Convert area into volume
AA = size(Hog1SignalingData,2);
for ii = 1:48;
    Hog1SignalingData(ii).Volm = (Hog1SignalingData2(ii).Volm).^1.5;
end;


%% figure handle
close; figure(11); set(gcf, 'Units', 'centimeters', 'Position', [0 0 21.5 25], 'PaperUnits', 'centimeters', 'PaperSize', [21.5 25]);
dx = 0.035; dy = 0.025; cmap = jet(7); lw=2; cmap([3 5],:)=.95*cmap([3 5],:); cmap = [cmap; hex2rgb('800080')]; cmap = [0,0,0;cmap]; 
  
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

%% 1, stimuli 
conc = [.1 .2 .4 .6 .8]; 
ID0=[2 3 5 7 9];
Hog1Conc.time = NaN(2501,6,8);
Hog1Conc.Conc = NaN(2501,6,8);
for j=1:5
    hh=subplotHJ(6,5,j,dy,dx); cla(hh); hold on; grid on;  box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
    count = 1; 
    data = DataSet1{100}; 
    plot(data.tt,-.1+data.stimulus_profile,'LineWidth', lw, 'color',cmap(count,:));
    count = count + 1; 
%     IDS  = [7:11:73]; 
    IDS  = [ID0(j):11:ID0(j)+11*6 ID0(j)+88];  
    for ii=1:size(IDS,2) 
        i = IDS(ii);
        data = DataSet1{i}; 
        timepoints = [data.tt(1):6/60:data.tt(end)]';% every 6 seconds (based on 25 minutes over 251 syringe pump phases)
        Hog1Conc.time(:,j,ii) = timepoints;
        Hog1Conc.Conc(:,j,ii) = -.1+data.stimulus(timepoints);
        plot(timepoints,-.1+data.stimulus(timepoints),'LineWidth', lw, 'color',cmap(count,:)); 
        count = count + 1; 
    end
   xlim([-3.5 50.5]); xticks([0:10:50]); ylim([.0 conc(j)]); yticks([0:.1:1]); 
end
N1 = numel(IDS); N2 = numel(ID0); % number of final concentrations and profiles
%%  2, volume data
Experiments = {'steps $ (t^{0}) $', 'root $ (\sqrt{t}) $', 'linears $ (t^{1}) $', 'quadratics $ (t^{2}) $', 'cubics $ (t^{3}) $', 'quintics $ (t^{5}) $', 'heptics $ (t^{7}) $', 'powers $ (a^t) $', 'sum polys'};
finalConcs = {[0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 1.0], [0.10 0.20:0.2:0.8], [[0.10 0.20:0.2:0.8]], ...
    [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]]}; 

concentrations = {'0.1M', '0.2M', '0.4M', '0.6M', '0.8M'}; 
ExperimentsID = {[1:3 6 8], [10:14], [15:19], [20:24], [25:29], [30:34], [35:39], [40:44]}; 

dd = 1;
Hog1Volume.Vol = NaN(1000,N1+1,N2); 
Hog1Volume.time = NaN(1000,N1+1,N2);

for i=1:N2%size(ExperimentsID,2)
    hh=subplotHJ(6,5,5+i,dy,dx); cla(hh); hold on; grid on; box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
%     text(25,1.25, concentrations{i}, 'Interpreter','latex', 'HorizontalAlignment', 'c');     

    % control 
    count = 1;     
    exp = 48; 
%     n_BRs = size(Hog1SignalingData(exp).Sel,2); 
    [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
    shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
    
    Hog1Volume.Vol(1:numel(dp),count,i) = Hog1SignalingData(exp).Volm(dp);
    Hog1Volume.time(1:numel(dp),count,i) = Hog1SignalingData(exp).tt(dp); 
    Hog1Volume.TimeDelay(count,i) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 
    count = count + 1; 

    for j=1:N1 %[1:4 5 6:7]
        exp = ExperimentsID{j}(i);
%         n_BRs = size(Hog1SignalingData(exp).Sel,2); 
        [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
        shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
        
        Hog1Volume.Vol(1:numel(dp),count,i) = Hog1SignalingData(exp).Volm(dp);
        Hog1Volume.time(1:numel(dp),count,i) = Hog1SignalingData(exp).tt(dp); 
        Hog1Volume.TimeDelay(count,i) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 
        count = count + 1; 
    end
    xlim([-3.5 50.5]); xticks([0:10:50]); % ylim([.72 1.27]); yticks([0:.2:1]); 
end
%% volume change 
time = Hog1Volume.time; 
time(:,3,4) = time(:,3,4)+.833333;  % shift

Vol0 = repmat(Hog1Volume.Vol(:,1,1),1,N1+1,N2)-Hog1Volume.Vol; % volume reduction
VolumeShrink = 1-Hog1Volume.Vol; 

%% 3, volume shrink 1, fit final change
    N_fits = 30; 
    best_params_1 = NaN(N_fits,N2,N1,2);  
    
    for i=1:N2
        Vol = Vol0(:,:,i); volShrink=VolumeShrink(:,:,i);
        hh=subplotHJ(6,5,10+i,dy,dx);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 25 25 25 25 25 25 25];
        for j=1:N1+1; hold on
           volShrinkPos=volShrink(:,j); 
           total_Vol_shrink(i,j) = sum(volShrinkPos(volShrinkPos>0)); 
            index = find(time(:,j,i)<=duration(j)+3); 
            this_data.time = time(index,j,i); this_data.vol = Vol(index,j);
            plot(time(index,j,i), Vol(index,j),'linewidth',1.5,'color',cmap(j,:)); 
            if j>=3
                f = 1;
                while f<=N_fits
                    try
                        [error, best_parameters, Vol_fit] = polynomial_function_fit(1,this_data,[poly_order(j) duration(j)]);
                        if f==N_fits
                        plot(time(index,j,i), Vol_fit,'linewidth',1.2,'color',cmap(j,:)); 
                        end
                        best_params_1(f,i,j-1,[1 2]) = best_parameters([1 2]);
                        f = f+1; 
                    catch
    %                     f = f+1;
                    end
                end
            end
        end
        box on; xlim([-3.5 27.5]); xticks([0:5:50]); %ylim([-.02 .32]); 
    end
    
%% 4, volume shrink 2, fit final change and polynomial order
    best_params_2 = NaN(N_fits,N2,N1,2);  
    
    for i=1:N2
        Vol = Vol0(:,:,i); 
        hh=subplotHJ(6,5,15+i,dy,dx);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 25 25 25 25 25 25 25];
        for j=1:N1+1; hold on
            index = find(time(:,j,i)<=duration(j)+3); 
            this_data.time = time(index,j,i); this_data.vol = Vol(index,j);
            plot(time(index,j,i), Vol(index,j),'linewidth',1.5,'color',cmap(j,:)); 
            f = 1;
            if j>=3
                while f<=N_fits
                    try
                        [error, best_parameters, Vol_fit] = polynomial_function_fit(3,this_data,[poly_order(j) duration(j)]);
                        if f==N_fits
                        plot(time(index,j,i), Vol_fit,'linewidth',1.2,'color',cmap(j,:)); 
                        end
                        best_params_2(f,i,j-1,[1 2]) = best_parameters([1 2]);
                        f = f+1; 
                    catch
    %                     f = f+1;
                    end
                end
            end
        end
        box on; xlim([-3.5 27.5]); xticks([0:5:50]); %ylim([-.02 .32]); 
    end 
    
 %% 5, volume shrink max value, data and estimated values
%  dx = 0.035; dy = 0.025;  
    for i=1:N2
        hh=subplotHJ(6,5,20+i,dy,dx);cla(hh); hold on; grid on; box on; 
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        Vol = Vol0(:,:,i); 
        
        index = find(time(:,1,i)<=25+3); 
        Vol25m = Vol(index,:); 
        VS = nanmax(Vol25m);
%         VS = VS./VS(2); 
        pars = reshape(best_params_2(:,i,:,1),N_fits,N1); 

        for j=1:size(VS,2)
            bar(j, VS(j), 'FaceColor',cmap(j,:),'EdgeColor','none');  
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=2
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)),'*','MarkerEdgeColor',.2*[1 1 1],'MarkerFaceColor',cmap(j,:),'MarkerSize',10); %cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
%         boxplot(reshape(best_params_1(:,i,:,1),N_fits,7),'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]);xticklabels([]); ylim([.0 .55]); yticks([0:.1:.5]); % 0:1:N1+1 ylim([.0 .33])
    end
    
%% 6, polynomial order, data and estimated values
    poly_order_vol = [];  
    for i=1:N2
        hh=subplotHJ(6,5,25+i,dy,dx); cla(hh); hold on; grid on; box on;         
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        
        pars = reshape(best_params_2(:,i,:,2),N_fits,N1);
        mu = round(nanmean(pars),2); 
        std = round(nanstd(pars),2)
        poly_order_vol = [poly_order_vol; mu];  
        for j=1:size(poly_order,2)
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=3
%                 bar(j, poly_order(j), 'FaceColor',cmap(j,:));  
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)),'*', 'color',cmap(j,:),'MarkerSize',10); %cmap(j,:));
%                 text(j,mu(j-1)*1.1,num2str(mu(j-1)),'color',cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact','OutlierSize',1)      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]); xticklabels([]); ylim([-.5 7.5]); yticks(poly_order(3:end-1));% yticks([0:1:30]); % 0:1:N1+1
    end
    Ma = nanmean(pars,1)    
    Sa = nanstd(pars,1)

    save('poly_order_vol_GN','poly_order_vol')
    save('total_Vol_shrink_GN','total_Vol_shrink')
%% save figs
% set(findall(gcf,'-property','FontSize'),'FontSize',7, 'color',.95*[1 1 1], 'defaultTextFontSize',7, {'DefaultAxesXColor','DefaultAxesYColor'},{.5*[1 1 1],.5*[1 1 1]}, 'FontName', 'Helvetica')
set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica')

figname = [dir_name, '/FigureS03GN'];
print(figname,'-depsc', '-r600');
%% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot individual figures:

%% 1, stimuli 
conc = [.1 .2 .4 .6 .8]; 
ID0=[2 3 5 7 9];
figure(2); clf;
for j=1:5
    hh=subplot(1,5,j); cla(hh); hold on; grid on;  box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
    count = 1; 
    data = DataSet1{100}; 
    plot(data.tt,-.1+data.stimulus_profile,'LineWidth', lw, 'color',cmap(count,:));
    count = count + 1; 
%     IDS  = [7:11:73]; 
    IDS  = [ID0(j):11:ID0(j)+11*6 ID0(j)+88];  
    for i=IDS 
        data = DataSet1{i}; 
        timepoints = [data.tt(1):6/60:data.tt(end)]';% every 6 seconds (based on 25 minutes over 251 syringe pump phases)
        plot(timepoints,-.1+data.stimulus(timepoints),'LineWidth', lw, 'color',cmap(count,:)); 
        count = count + 1; 
    end
   xlim([-3.5 50.5]); xticks([0:10:50]); ylim([.0 conc(j)]); yticks([0:.1:1]); 
end
N1 = numel(IDS); N2 = numel(ID0); % number of final concentrations and profiles

%%  2, volume data
Experiments = {'steps $ (t^{0}) $', 'root $ (\sqrt{t}) $', 'linears $ (t^{1}) $', 'quadratics $ (t^{2}) $', 'cubics $ (t^{3}) $', 'quintics $ (t^{5}) $', 'heptics $ (t^{7}) $', 'powers $ (a^t) $', 'sum polys'};
finalConcs = {[0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 1.0], [0.10 0.20:0.2:0.8], [[0.10 0.20:0.2:0.8]], ...
    [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]]}; 

concentrations = {'0.1M', '0.2M', '0.4M', '0.6M', '0.8M'}; 
ExperimentsID = {[1:3 6 8], [10:14], [15:19], [20:24], [25:29], [30:34], [35:39], [40:44]}; 

dd = 1;
Hog1Volume.Vol = NaN(1000,N1+1,N2); 
Hog1Volume.time = NaN(1000,N1+1,N2); 

figure(3); clf;
for i=1:N2%size(ExperimentsID,2)
    hh=subplot(1,5,i); cla(hh); hold on; grid on; box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
%     text(25,1.25, concentrations{i}, 'Interpreter','latex', 'HorizontalAlignment', 'c');     

    % control 
    count = 1;     
    exp = 48; 
%     n_BRs = size(Hog1SignalingData(exp).Sel,2); 
    [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
    shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
    
    Hog1Volume.Vol(1:numel(dp),count,i) = Hog1SignalingData(exp).Volm(dp);
    Hog1Volume.time(1:numel(dp),count,i) = Hog1SignalingData(exp).tt(dp); 
    Hog1Volume.TimeDelay(count,i) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 
    count = count + 1; 

    for j=1:N1 %[1:4 5 6:7]
        exp = ExperimentsID{j}(i);
%         n_BRs = size(Hog1SignalingData(exp).Sel,2); 
        [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
        shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
        
        Hog1Volume.Vol(1:numel(dp),count,i) = Hog1SignalingData(exp).Volm(dp);
        Hog1Volume.time(1:numel(dp),count,i) = Hog1SignalingData(exp).tt(dp); 
        Hog1Volume.TimeDelay(count,i) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 
        count = count + 1; 
    end
    xlim([-3.5 50.5]); xticks([0:10:50]); % ylim([.72 1.27]); yticks([0:.2:1]); 
end
%% volume change 
time = Hog1Volume.time; 
time(:,3,4) = time(:,3,4)+.833333;  % shift

Vol0 = repmat(Hog1Volume.Vol(:,1,1),1,N1+1,N2)-Hog1Volume.Vol; % volume reduction
VolumeShrink = 1-Hog1Volume.Vol; 

%% 3, volume shrink 1, fit final change
    N_fits = 30; 
    best_params_1 = NaN(N_fits,N2,N1,2);  
    
figure(4); clf;
    for i=1:N2
        Vol = Vol0(:,:,i); volShrink=VolumeShrink(:,:,i);
        hh=subplot(1,5,i);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 25 25 25 25 25 25 25];
        for j=1:N1+1; hold on
           volShrinkPos=volShrink(:,j); 
           total_Vol_shrink(i,j) = sum(volShrinkPos(volShrinkPos>0)); 
            index = find(time(:,j,i)<=duration(j)+3); 
            this_data.time = time(index,j,i); this_data.vol = Vol(index,j);
            plot(time(index,j,i), Vol(index,j),'linewidth',1.5,'color',cmap(j,:)); 
            if j>=3
                f = 1;
                while f<=N_fits
                    try
                        [error, best_parameters, Vol_fit] = polynomial_function_fit(1,this_data,[poly_order(j) duration(j)]);
                        if f==N_fits
                        plot(time(index,j,i), Vol_fit,'linewidth',1.2,'color',cmap(j,:)); 
                        end
                        best_params_1(f,i,j-1,[1 2]) = best_parameters([1 2]);
                        f = f+1; 
                    catch
    %                     f = f+1;
                    end
                end
            end
        end
        box on; xlim([-3.5 27.5]); xticks([0:5:50]); %ylim([-.02 .32]); 
    end
    
%% 4, volume shrink 2, fit final change and polynomial order
    best_params_2 = NaN(N_fits,N2,N1,2);  
    
figure(5); clf;
    for i=1:N2
        Vol = Vol0(:,:,i); 
        hh=subplot(1,5,i);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 25 25 25 25 25 25 25];
        for j=1:N1+1; hold on
            index = find(time(:,j,i)<=duration(j)+3); 
            this_data.time = time(index,j,i); this_data.vol = Vol(index,j);
            plot(time(index,j,i), Vol(index,j),'linewidth',1.5,'color',cmap(j,:)); 
            f = 1;
            if j>=3
                while f<=N_fits
                    try
                        [error, best_parameters, Vol_fit] = polynomial_function_fit(3,this_data,[poly_order(j) duration(j)]);
                        if f==N_fits
                        plot(time(index,j,i), Vol_fit,'linewidth',1.2,'color',cmap(j,:)); 
                        end
                        best_params_2(f,i,j-1,[1 2]) = best_parameters([1 2]);
                        f = f+1; 
                    catch
    %                     f = f+1;
                    end
                end
            end
        end
        box on; xlim([-3.5 27.5]); xticks([0:5:50]); %ylim([-.02 .32]); 
    end 
    
 %% 5, volume shrink max value, data and estimated values
%  dx = 0.035; dy = 0.025;  

figure(6); clf;
    for i=1:N2
        hh=subplot(1,5,i);cla(hh); hold on; grid on; box on; 
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        Vol = Vol0(:,:,i); 
        
        index = find(time(:,1,i)<=25+3); 
        Vol25m = Vol(index,:); 
        VS = nanmax(Vol25m);
%         VS = VS./VS(2); 
        pars = reshape(best_params_2(:,i,:,1),N_fits,N1); 

        for j=1:size(VS,2)
            bar(j, VS(j), 'FaceColor',cmap(j,:),'EdgeColor','none');  
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=2
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)),'*','MarkerEdgeColor',.2*[1 1 1],'MarkerFaceColor',cmap(j,:),'MarkerSize',10); %cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
%         boxplot(reshape(best_params_1(:,i,:,1),N_fits,7),'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]);xticklabels([]); ylim([.0 .33]); yticks([0:.1:.30]); % 0:1:N1+1
    end
    
%% 6, polynomial order, data and estimated values
    poly_order_vol = [];

figure(7); clf;
    for i=1:N2
        hh=subplot(1,5,i); cla(hh); hold on; grid on; box on;         
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        
        pars = reshape(best_params_2(:,i,:,2),N_fits,N1);
        mu = round(nanmean(pars),2); poly_order_vol = [poly_order_vol; mu];  
        for j=1:size(poly_order,2)
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=3
%                 bar(j, poly_order(j), 'FaceColor',cmap(j,:));  
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)),'*', 'color',cmap(j,:),'MarkerSize',10); %cmap(j,:));
%                 text(j,mu(j-1)*1.1,num2str(mu(j-1)),'color',cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact','OutlierSize',1)      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]); xticklabels([]); ylim([-.5 7.5]); yticks(poly_order(3:end-1));% yticks([0:1:30]); % 0:1:N1+1
    end
Mv = nanmean(pars,1)    
Sv = nanstd(pars,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot individual figures for area converted to volume

%% Convert area into volume
AA = size(Hog1SignalingData,2);
for ii = 1:48;
    Hog1SignalingData(ii).Volm = (Hog1SignalingData2(ii).Volm).^1.5;
end;

%% 1, stimuli 
conc = [.1 .2 .4 .6 .8]; 
ID0=[2 3 5 7 9];
figure(12); clf;
for j=1:5
    hh=subplot(1,5,j); cla(hh); hold on; grid on;  box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
    count = 1; 
    data = DataSet1{100}; 
    plot(data.tt,-.1+data.stimulus_profile,'LineWidth', lw, 'color',cmap(count,:));
    count = count + 1; 
%     IDS  = [7:11:73]; 
    IDS  = [ID0(j):11:ID0(j)+11*6 ID0(j)+88];  
    for i=IDS 
        data = DataSet1{i}; 
        timepoints = [data.tt(1):6/60:data.tt(end)]';% every 6 seconds (based on 25 minutes over 251 syringe pump phases)
        plot(timepoints,-.1+data.stimulus(timepoints),'LineWidth', lw, 'color',cmap(count,:)); 
        count = count + 1; 
    end
   xlim([-3.5 50.5]); xticks([0:10:50]); ylim([.0 conc(j)]); yticks([0:.1:1]); 
end
N1 = numel(IDS); N2 = numel(ID0); % number of final concentrations and profiles

%%  2, volume data
Experiments = {'steps $ (t^{0}) $', 'root $ (\sqrt{t}) $', 'linears $ (t^{1}) $', 'quadratics $ (t^{2}) $', 'cubics $ (t^{3}) $', 'quintics $ (t^{5}) $', 'heptics $ (t^{7}) $', 'powers $ (a^t) $', 'sum polys'};
finalConcs = {[0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 1.0], [0.10 0.20:0.2:0.8], [[0.10 0.20:0.2:0.8]], ...
    [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]], [[0.10 0.20:0.2:0.8]]}; 

concentrations = {'0.1M', '0.2M', '0.4M', '0.6M', '0.8M'}; 
ExperimentsID = {[1:3 6 8], [10:14], [15:19], [20:24], [25:29], [30:34], [35:39], [40:44]}; 

dd = 1;
Hog1Volume.Vol = NaN(1000,N1+1,N2); 
Hog1Volume.time = NaN(1000,N1+1,N2); 

figure(13); clf;
for i=1:N2%size(ExperimentsID,2)
    hh=subplot(1,5,i); cla(hh); hold on; grid on; box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
%     text(25,1.25, concentrations{i}, 'Interpreter','latex', 'HorizontalAlignment', 'c');     

    % control 
    count = 1;     
    exp = 48; 
%     n_BRs = size(Hog1SignalingData(exp).Sel,2); 
    [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
    shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
    
    Hog1Volume.Vol(1:numel(dp),count,i) = Hog1SignalingData(exp).Volm(dp);
    Hog1Volume.time(1:numel(dp),count,i) = Hog1SignalingData(exp).tt(dp); 
    Hog1Volume.TimeDelay(count,i) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 
    count = count + 1; 

    for j=1:N1 %[1:4 5 6:7]
        exp = ExperimentsID{j}(i);
%         n_BRs = size(Hog1SignalingData(exp).Sel,2); 
        [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData);          
        shadedErrorBar(Hog1SignalingData(exp).tt(dp),Hog1SignalingData(exp).Volm(dp),Hog1SignalingData(exp).Vols(dp),{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
        
        Hog1Volume.Vol(1:numel(dp),count,i) = Hog1SignalingData(exp).Volm(dp);
        Hog1Volume.time(1:numel(dp),count,i) = Hog1SignalingData(exp).tt(dp); 
        Hog1Volume.TimeDelay(count,i) = numel(find(Hog1SignalingData(exp).tt(dp)<=0)); 
        count = count + 1; 
    end
    xlim([-3.5 50.5]); xticks([0:10:50]); % ylim([.72 1.27]); yticks([0:.2:1]); 
end
%% volume change 
time = Hog1Volume.time; 
time(:,3,4) = time(:,3,4)+.833333;  % shift

Vol0 = repmat(Hog1Volume.Vol(:,1,1),1,N1+1,N2)-Hog1Volume.Vol; % volume reduction
VolumeShrink = 1-Hog1Volume.Vol; 

%% 3, volume shrink 1, fit final change
    N_fits = 30; 
    best_params_1 = NaN(N_fits,N2,N1,2);  
    
figure(14); clf;
    for i=1:N2
        Vol = Vol0(:,:,i); volShrink=VolumeShrink(:,:,i);
        hh=subplot(1,5,i);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 25 25 25 25 25 25 25];
        for j=1:N1+1; hold on
           volShrinkPos=volShrink(:,j); 
           total_Vol_shrink(i,j) = sum(volShrinkPos(volShrinkPos>0)); 
            index = find(time(:,j,i)<=duration(j)+3); 
            this_data.time = time(index,j,i); this_data.vol = Vol(index,j);
            plot(time(index,j,i), Vol(index,j),'linewidth',1.5,'color',cmap(j,:)); 
            if j>=3
                f = 1;
                while f<=N_fits
                    try
                        [error, best_parameters, Vol_fit] = polynomial_function_fit(1,this_data,[poly_order(j) duration(j)]);
                        if f==N_fits
                        plot(time(index,j,i), Vol_fit,'linewidth',1.2,'color',cmap(j,:)); 
                        end
                        best_params_1(f,i,j-1,[1 2]) = best_parameters([1 2]);
                        f = f+1; 
                    catch
    %                     f = f+1;
                    end
                end
            end
        end
        box on; xlim([-3.5 27.5]); xticks([0:5:50]); %ylim([-.02 .32]); 
    end
    
%% 4, volume shrink 2, fit final change and polynomial order
    best_params_2 = NaN(N_fits,N2,N1,2);  
    
figure(15); clf;
    for i=1:N2
        Vol = Vol0(:,:,i); 
        hh=subplot(1,5,i);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 25 25 25 25 25 25 25];
        for j=1:N1+1; hold on
            index = find(time(:,j,i)<=duration(j)+3); 
            this_data.time = time(index,j,i); this_data.vol = Vol(index,j);
            plot(time(index,j,i), Vol(index,j),'linewidth',1.5,'color',cmap(j,:)); 
            f = 1;
            if j>=3
                while f<=N_fits
                    try
                        [error, best_parameters, Vol_fit] = polynomial_function_fit(3,this_data,[poly_order(j) duration(j)]);
                        if f==N_fits
                        plot(time(index,j,i), Vol_fit,'linewidth',1.2,'color',cmap(j,:)); 
                        end
                        best_params_2(f,i,j-1,[1 2]) = best_parameters([1 2]);
                        f = f+1; 
                    catch
    %                     f = f+1;
                    end
                end
            end
        end
        box on; xlim([-3.5 27.5]); xticks([0:5:50]); %ylim([-.02 .32]); 
    end 
    
 %% 5, volume shrink max value, data and estimated values
%  dx = 0.035; dy = 0.025;  

figure(16); clf;
    for i=1:N2
        hh=subplot(1,5,i);cla(hh); hold on; grid on; box on; 
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        Vol = Vol0(:,:,i); 
        
        index = find(time(:,1,i)<=25+3); 
        Vol25m = Vol(index,:); 
        VS = nanmax(Vol25m);
%         VS = VS./VS(2); 
        pars = reshape(best_params_2(:,i,:,1),N_fits,N1); 

        for j=1:size(VS,2)
            bar(j, VS(j), 'FaceColor',cmap(j,:),'EdgeColor','none');  
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=2
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)),'*','MarkerEdgeColor',.2*[1 1 1],'MarkerFaceColor',cmap(j,:),'MarkerSize',10); %cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
%         boxplot(reshape(best_params_1(:,i,:,1),N_fits,7),'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]);xticklabels([]); ylim([.0 .33]); yticks([0:.1:.30]); % 0:1:N1+1
    end
    
%% 6, polynomial order, data and estimated values
    poly_order_vol = [];

figure(17); clf;
    for i=1:N2
        hh=subplot(1,5,i); cla(hh); hold on; grid on; box on;         
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        
        pars = reshape(best_params_2(:,i,:,2),N_fits,N1);
        mu = round(nanmean(pars),2); poly_order_vol = [poly_order_vol; mu];  
        for j=1:size(poly_order,2)
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=3
%                 bar(j, poly_order(j), 'FaceColor',cmap(j,:));  
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)),'*', 'color',cmap(j,:),'MarkerSize',10); %cmap(j,:));
%                 text(j,mu(j-1)*1.1,num2str(mu(j-1)),'color',cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact','OutlierSize',1)      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]); xticklabels([]); ylim([-.5 7.5]); yticks(poly_order(3:end-1));% yticks([0:1:30]); % 0:1:N1+1
    end

%% Compute AUC
%% 1, stimuli 
conc = [.1 .2 .4 .6 .8];
poly = {'step','x^1^/^2','x^1','x^2','x^3','x^4','x^7','exp'}; 
ID0=[2 3 5 7 9];
figure(12); clf;
j = 1;
ii = 1;
for j=1:5
    hh=subplot(1,5,j); cla(hh); hold on; grid on;  box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
    count = 1; 
    data = DataSet1{100}; 
    timepoints = [data.tt(1):6/60:data.tt(end)]';% every 6 seconds (based on 25 minutes over 251 syringe pump phases)
    timepoints1 = timepoints(1001:1501);
    CS = cumsum(-.1+data.stimulus(timepoints1));
    CSall(:,i) = CS;
    CSmax(j,1) = CS(501);
    ProfileConc{j,1} = {conc(j),poly{ii}};
    plot(timepoints1,CS,'LineWidth', lw, 'color',cmap(count,:)); 
    count = count + 1; 
%     IDS  = [7:11:73]; 
    IDS  = [ID0(j):11:ID0(j)+11*6 ID0(j)+88];  
    for ii=1:size(IDS,2);
        i = IDS(ii);
        data = DataSet1{i}; 
        timepoints = [data.tt(1):6/60:data.tt(end)]';% every 6 seconds (based on 25 minutes over 251 syringe pump phases)
        timepoints1 = timepoints(1001:1501);
        CS = cumsum(-.1+data.stimulus(timepoints1));
        CSall(:,j,ii) = CS;
        CSmax(j,ii) = CS(501);
        ProfileConc{j,ii} = [num2str(conc(j)) '-' poly{ii}];
        ProfileConc2{j,ii} = {conc(j),poly{ii}};
        plot(timepoints1,CS,'LineWidth', lw, 'color',cmap(count,:)); 
        count = count + 1; 
    end
    legend(ProfileConc{j,:},'Location','northwest');
    xlim([-3.5 50.5]); ylim([.0 max(CSall(501,IDS))]);xticks([0:10:50]); yticks([0:max(CSall(501,IDS))/10:max(CSall(501,IDS))]);  % ylim([.0 conc(j)]);
end
N1 = numel(IDS); N2 = numel(ID0); % number of final concentrations and profiles

%% compare volume change for similar total AUC of NaCl:
i = 1;
rowAll = NaN(15,5);
colAll = NaN(15,5);
figure(21); clf;
for i = 1:5;
    i
    AA = CSmax(i,1);
    [row, col] = find(CSmax >= 0.8*CSmax(i,1) & CSmax <= 1.2*CSmax(i,1));
    rowAll(1:size(row,1),i) = row;
    colAll(1:size(row,1),i) = col;
    subplot(2,5,i);
    count = 1;
    for j = 1:size(row,1);
        j
        subplot(2,5,i);
        ttt = Hog1Conc.time(:,row(j),col(j));
        Conc = Hog1Conc.Conc(:,row(j),col(j));
        plot(ttt,Conc,'LineWidth', lw); hold on; % , 'color',cmap(count,:)
        xlim([-3.5 50.5]);
        title(['AUC =' num2str(AA) '+/- 20%']);
        subplot(2,5,i+5);
        tv = Hog1Volume.time(:,col(j),row(j));
        Voll = Hog1Volume.Vol(:,col(j),row(j));
        plot(tv,Voll,'LineWidth', lw); hold on; % , 'color',cmap(count,:)
        xlim([-3.5 50.5]);
        count = count + 1;

    end;
    legend(ProfileConc{row,col},'Location','northwest');
end;
    

%% Need to update for Hog1 signaling
%% compare Hog1 change for similar total AUC of NaCl:
i = 1;
rowAll = NaN(15,5);
colAll = NaN(15,5);
figure(21); clf;
for i = 1:5;
    i
    AA = CSmax(i,1);
    [row, col] = find(CSmax >= 0.8*CSmax(i,1) & CSmax <= 1.2*CSmax(i,1));
    rowAll(1:size(row,1),i) = row;
    colAll(1:size(row,1),i) = col;
    subplot(2,5,i);
    count = 1;
    for j = 1:size(row,1);
        j
        % NaCl profiles
        subplot(4,5,i);
        ttt = Hog1Conc.time(:,row(j),col(j));
        Conc = Hog1Conc.Conc(:,row(j),col(j));
        plot(ttt,Conc,'LineWidth', lw); hold on; % , 'color',cmap(count,:)
        xlim([-3.5 50.5]);
        title(['AUC =' num2str(AA) '+/- 20%']);
        
        % Integrated NaCl = AUC
        subplot(4,5,i+10);
        
        % Volume profiles
        subplot(4,5,i+10);
        tv = Hog1Volume.time(:,col(j),row(j));
        Voll = Hog1Volume.Vol(:,col(j),row(j));
        plot(tv,Voll,'LineWidth', lw); hold on; % , 'color',cmap(count,:)
        xlim([-3.5 50.5]);
        
        % Hog1 profiles
        subplot(4,5,i+5);
        tv = Hog1Volume.time(:,col(j),row(j));
        Voll = Hog1Volume.Vol(:,col(j),row(j));
        plot(tv,Voll,'LineWidth', lw); hold on; % , 'color',cmap(count,:)
        xlim([-3.5 50.5]);
        count = count + 1;
    end;
    subplot(2,5,i)
    legend(ProfileConc{row,col},'Location','northwest');
end;
