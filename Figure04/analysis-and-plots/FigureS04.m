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

%% figure handle
close; figure(2); set(gcf, 'Units', 'centimeters', 'Position', [0 0 21.5 25], 'PaperUnits', 'centimeters', 'PaperSize', [21.5 25]);
dx = 0.035; dy = 0.025; cmap = jet(7); lw=2; cmap([3 5],:)=.95*cmap([3 5],:); cmap = [cmap; hex2rgb('800080')]; cmap = [0,0,0;cmap]; 
  
set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

%% 1, stimuli 
conc = [.1 .2 .4 .6 .8]; 
ID0=[2 3 5 7 9];
for j=1:5
    hh=subplotHJ(6,5,j,dy,dx); cla(hh); hold on; grid on;  box on; 
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
Hog1Volume.Hog1 = NaN(1000,N1+1,N2); 
Hog1Volume.time = NaN(1000,N1+1,N2); 

for i=1:N2%size(ExperimentsID,2)
    hh=subplotHJ(6,5,5+i,dy,dx); cla(hh); hold on; grid on; box on; 
    set(groot,'defaultAxesTickLabelInterpreter','tex');
%     text(25,1.25, concentrations{i}, 'Interpreter','latex', 'HorizontalAlignment', 'c');     

    % control 
    count = 1;   
    data = DataSet1{100}; 
    shadedErrorBar(data.tt,data.mHog1nuc,data.sHog1nuc,{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
    Hog1Volume.Hog1(1:numel(data.mHog1nuc),count,i) = data.mHog1nuc;
    Hog1Volume.time(1:numel(data.tt),count,i) = data.tt;     
    count = count + 1; 
%     IDS  = [7:11:73]; 
    IDS  = [ID0(i):11:ID0(i)+11*6 ID0(i)+88];  
    
    for j=IDS %[1:4 5 6:7]        
        
        data = DataSet1{j}; 
        shadedErrorBar(data.tt,data.mHog1nuc,data.sHog1nuc,{'LineWidth', lw, 'color',cmap(count,:)}, 0.2); 
        Hog1Volume.Hog1(1:numel(data.mHog1nuc),count,i) = data.mHog1nuc;
        Hog1Volume.time(1:numel(data.tt),count,i) = data.tt;            
        count = count + 1; 
    end
    xlim([-3.5 50.5]); xticks([0:10:50]); % ylim([.72 1.27]); yticks([0:.2:1]); 
end
%% volume change and Hog1 versus volume
time = Hog1Volume.time; 
% time(:,3,4) = time(:,3,4)+.833333;  % shift

Hog1 = Hog1Volume.Hog1; 

%% 3, volume shrink 1, fit final change
    N_fits = 30; 
    best_params_1 = NaN(N_fits,N2,N1,2);  
    
    for i=1:N2
        Hog1_ = Hog1(:,:,i); 
        hh=subplotHJ(6,5,10+i,dy,dx);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 6 15 25 25 25 25 15];
        for j=1:N1+1; hold on
            index = find(time(:,j,i)>=-5 & time(:,j,i)<=duration(j)+2); 
            this_data.time = time(index,j,i); this_data.vol = Hog1_(index,j);
            if j>=3
                plot(time(index,j,i), Hog1_(index,j),'linewidth',1.5,'color',cmap(j,:)); 
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
        Hog1_ = Hog1(:,:,i); 
        hh=subplotHJ(6,5,15+i,dy,dx);cla(hh); hold on; grid on
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        poly_order = [NaN .01 .5 1 2 3 5 7 1];
        duration = [NaN 1 6 15 25 25 25 25 15];
        for j=1:N1+1; hold on
            index = find(time(:,j,i)>=-5 & time(:,j,i)<=duration(j)+2); 
            this_data.time = time(index,j,i); this_data.vol = Hog1_(index,j);
            f = 1;
            if j>=3
                plot(time(index,j,i), Hog1_(index,j),'linewidth',1.5,'color',cmap(j,:)); 
                while f<=N_fits
                    try
                        [error, best_parameters, Vol_fit] = polynomial_function_fit(3,this_data,[poly_order(j) duration(j)]);
                        if f==N_fits
                        plot(time(index,j), Vol_fit,'linewidth',1.2,'color',cmap(j,:)); 
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
        Hog1_ = Hog1(:,:,i); 
        index = find(time(:,1,i)>=-5 & time(:,1,i)<=25+3); 
        Vol25m = Hog1_(index,:); 
        VS = nanmax(Vol25m);
%         VS = VS./VS(2); 
        pars = reshape(best_params_2(:,i,:,1),N_fits,N1); 

        for j=1:size(VS,2)
            bar(j, VS(j), 'FaceColor',cmap(j,:),'EdgeColor','none');  
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=2
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)), '*','MarkerEdgeColor',.2*[1 1 1],'MarkerFaceColor',cmap(j,:),'MarkerSize',10); %cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
%         boxplot(reshape(best_params_1(:,i,:,1),N_fits,7),'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact')      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]);xticklabels([]); ylim([.0 1.5]); yticks([0:.5:1.50]); % 0:1:N1+1
    end
    
%% 6, polynomial order, data and estimated values
    poly_order_Hog1 = [];  

    for i=1:N2
        hh=subplotHJ(6,5,25+i,dy,dx); cla(hh); hold on; grid on; box on;         
        set(groot,'defaultAxesTickLabelInterpreter','tex');
        
        pars = reshape(best_params_2(:,i,:,2),N_fits,N1);
        mu = round(nanmean(pars),2); poly_order_Hog1 = [poly_order_Hog1; mu];
        
        for j=1:size(poly_order,2)
%             plot(j, VS(j), '*', 'color',cmap(j,:)); 
            if j>=3
%                 bar(j, poly_order(j), 'FaceColor',cmap(j,:));  
                errorbar(j, nanmean(pars(:,j-1)), nanstd(pars(:,j-1)),'color', .2*[1 1 1]); %cmap(j,:)); 
                plot(j, nanmean(pars(:,j-1)), '*','color',cmap(j,:),'MarkerSize',10); %cmap(j,:));
            end
        end
%         boxplot(pars,'Positions',[2:8],'Colors',cmap(2:end,:),'PlotStyle', 'compact','OutlierSize',1)      
        xlim([.5 N1+1.5]); xticks([0:1:N1+1]); xticklabels([]); ylim([-.5 10.5]); yticks(poly_order(3:end-1));% yticks([0:1:30]); % 0:1:N1+1
    end
        save('poly_order_Hog1','poly_order_Hog1')

%% save figs
% set(findall(gcf,'-property','FontSize'),'FontSize',7, 'color',.95*[1 1 1], 'defaultTextFontSize',7, {'DefaultAxesXColor','DefaultAxesYColor'},{.5*[1 1 1],.5*[1 1 1]}, 'FontName', 'Helvetica')
set(findall(gcf,'-property','FontSize'),'FontSize',7, 'defaultTextFontSize',7, 'FontName', 'Helvetica')

figname = [dir_name, '/FigureS04'];
print(figname,'-depsc', '-r600');
%% end
