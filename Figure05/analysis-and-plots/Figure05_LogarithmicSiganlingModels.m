function LogarithmicSiganlingData =  Figure05_LogarithmicSiganlingModels() 
clearvars -except Hog1SignalingData; 

dir_name = ['plots'];
mkdir(dir_name); 

col_csle(1,:) = .6*[1 1 1]; %hex2rgb('f0f5f9'); %ctrl
col_csle(2,:) = hex2rgb('090707'); %step black
% col_csle(2,:) = hex2rgb('87314e'); %step red
col_csle(3,:) = hex2rgb('00b9f1'); %linear
col_csle(4,:) = hex2rgb('A593E0'); %linear
col_csle(5,:) = hex2rgb('87314e'); % sum
col_csle(6,:) = [0 0 1]; 
%% define stimuli profiles 

%control stimuli u=u0;
constant_stimuli_input = @(t, conc_start) (t-t + conc_start);

% stimuli profile as polynomial functions
polynomial_stimuli_profile = @(t, conc_start, k, duration, conc_change) max(conc_start,min([conc_start + (conc_change.*t.^k/(duration)^k),conc_start+conc_change]));
% polynomial_stimuli_profile = @(t, conc_start, k, duration, conc_change) conc_start + (conc_change.*t.^k/(duration)^k);

% stimuli profile as exponential functions
exponential_stimuli_profile = @(t, conc_start, duration, conc_change) max(conc_start, min([(conc_start -1 + ((1+conc_change)^(1/(duration)))^t),conc_start+conc_change]));

% stairs stimuli
sympref('HeavisideAtOrigin',1);
stair_stimuli_profile = @(t, conc_start, duration, conc_change) conc_start + conc_change * heaviside(t - .001) ...
    + conc_change * heaviside(t - duration - .001) ...
    + conc_change * heaviside(t - 2*duration - .001) ... 
    + conc_change * heaviside(t - 3*duration - .001) ... 
    + conc_change * heaviside(t - 4*duration - .001);

% stairs stimuli
sympref('HeavisideAtOrigin',1);  ramp = 5; 
stair_stimuli_profile2 = @(t, conc_start, duration, conc_change) conc_start + polynomial_stimuli_profile(t, 0, 1, ramp, conc_change) * heaviside(t - .001) ...
    + polynomial_stimuli_profile(t-duration, 0, 1, ramp, .65*conc_change) * heaviside(t - duration - .001) ...
    + polynomial_stimuli_profile(t-2*duration, 0, 1, ramp, .35*conc_change) * heaviside(t - 2*duration - .001);

% pulse stimuli 
pulsatile_stimuli_profile = @(t, conc_start, duration, conc_change) conc_start + .5*conc_change*(square((2*pi/duration)*(t-1))+1); 

% stimuli profiles; constant, stairs, linear, exponential
start_concnetration = 1; 
stimuli_duration = 50; 
stair_step_changes = 1/3; 
pulse_step_changes = 1/3; 
conc_change = 1; 
w = 5; 

control_stimuli = @(t) constant_stimuli_input(t, start_concnetration); 
stair_stimuli = @(t) stair_stimuli_profile(t, start_concnetration, stimuli_duration, stair_step_changes); 
linear_stimuli = @(t) polynomial_stimuli_profile(t, start_concnetration, 1, stimuli_duration, conc_change); 
exponential_stimuli = @(t) exponential_stimuli_profile(t, start_concnetration, stimuli_duration, conc_change); 
pulsatile_stimuli = @(t) pulsatile_stimuli_profile(t, start_concnetration, stimuli_duration, pulse_step_changes); 

a = (1+conc_change)^(1/stimuli_duration);
poly_stimuli = @(t,k) ((stimuli_duration*log(a)).^k)*(conc_change/(stimuli_duration^k))*(1/factorial(k))*((stimuli_duration^k)/conc_change)*(polynomial_stimuli_profile(t, start_concnetration, k, stimuli_duration, conc_change)-start_concnetration); 

stimuli_input{1} = @(t) constant_stimuli_input(t, start_concnetration); 
stimuli_input{2} = @(t) stair_stimuli_profile(t, start_concnetration, 20, stair_step_changes); 
stimuli_input{3} = @(t) polynomial_stimuli_profile(t, start_concnetration, 1, stimuli_duration, conc_change); 
stimuli_input{4} = @(t) exponential_stimuli_profile(t, start_concnetration, stimuli_duration, conc_change); 
stimuli_input{5} = @(t) poly_stimuli(t,1)+poly_stimuli(t,2)+poly_stimuli(t,3)+poly_stimuli(t,5)+poly_stimuli(t,7)+start_concnetration;
stimuli_input{6} = @(t) pulsatile_stimuli_profile(t, start_concnetration, 20, pulse_step_changes); 

%  simulate LOGARITHMIC SIGNALING responses
close; figure();clf; set(gcf, 'Units', 'centimeters', 'Position', [0 0 18 9], 'PaperUnits', 'centimeters', 'PaperSize', [18 9]); 
dx = 0.008; dy = 0.015; lw=2; 

hold on
timepoints4 = {[-10:.1:60]', [-10:.5:60]', [-10:1.5:60]', [-10:1.5:60]',[-10:1.5:60]',[-10:.5:60]'}; 
timepoints =  timepoints4{1}; 
i=1;
for t=1:length(timepoints)
    stair(i) = stimuli_input{2}(timepoints(t));
    pulse(i) = stimuli_input{6}(timepoints(t));
    i=i+1;
end
coef = [.55 1 1 1 .55];  
count = 1; 
for j=[3 4 6 2 5]
    %% stimuli 
    timepoints =  timepoints4{j}; 
    N_TP = length(timepoints); 
    i=1; 

    for t=1:length(timepoints)
        stimuli(i,1) = stimuli_input{j}(timepoints(t)); 
        i=i+1;
    end
    stimuli = stimuli; 
    stimuli_rate = max(0,4*[0; diff(stimuli)]); 
    
    %% log signaling
    response0 = stimuli_rate./([1; stimuli(1:end-1)]); 

    for i=1:length(response0)
        responseLOG(i,1) = sum(response0([max(1,i-w):i])); 
    end
    
%     i=1; 
%     for t=1:length(timepoints)
%         stimuli(i,1) = stimuli_input{j}(timepoints(t)); 
%         i=i+1;
%     end
    
    LogDataSet{j}.stimuli_profile = stimuli_input{j}; 
    LogDataSet{j}.tt = timepoints; 
    LogDataSet{j}.stimuli = stimuli; 
    LogDataSet{j}.responseLOG = responseLOG; 
%     LogDataSet{j}.responseLIN = responseLIN; 

    subplotHJ(4,5,count,dy,dx); hold on; % stimuli 
    if j==2
        plot(timepoints4{1},stair, 'LineWidth', lw, 'color', col_csle(j,:)); 
    else
        plot(timepoints,stimuli, 'LineWidth', lw, 'color', col_csle(j,:)); 
    end
    xlim([-5 57]); xticks([0 50]);  
    ylim([.95 2.05]); yticks([1 2]); box on; 
    xticklabels([]); yticklabels([]); 

    subplotHJ(4,5,5*2+count,dy,dx); hold on; % logrithmic signaling
%     R_LOG = min(1.3,responseLOG.^1.7); 
    R_LOG = coef(j-1)*responseLOG.^1.7; 
    plot(timepoints,R_LOG, 'LineWidth', 2*lw, 'color', col_csle(j,:)); 
%     shadedErrorBar(timepoints,response,STD,{'LineWidth', 2}, .15); 
    xlim([-5 57]); xticks([0 50]);  
    ylim([-.03 1.03]); yticks([0 1]); box on
    xticklabels([]); yticklabels([]); 

    %% linear signaling
    response0 = stimuli_rate; 4*[0; diff(stimuli)];

    for i=1:length(response0)
        responseLIN(i,1) = sum(response0([max(1,i-w):i])); 
    end
    
    subplotHJ(4,5,5*1+count,dy,dx); hold on; % linear signaling
%     R_LIN = min(1.3,responseLIN.^1.7); 
    R_LIN = coef(j-1)*responseLIN.^1.7; 
    plot(timepoints,R_LIN, 'LineWidth', 2*lw, 'color', col_csle(j,:)); 
    xlim([-5 57]); xticks([0 50]);  
    ylim([-.03 1.03]); yticks([0 1]); box on
    xticklabels([]); yticklabels([]); 
    
    %% log + lin 
    R_LINLOG = .5*(R_LIN + R_LOG); 
    subplotHJ(4,5,5*3+count,dy,dx); hold on; % hybrid signaling
    plot(timepoints,R_LINLOG, 'LineWidth', 2*lw, 'color', col_csle(j,:)); 
    xlim([-5 57]); xticks([0 50]);  
    ylim([-.03 1.03]); yticks([0 1]); box on
    xticklabels([]); yticklabels([]); 
    
    clear stimuli responseLOG responseLIN 
    count = count +1; 
end

save([dir_name,'/LogDataSet'],'LogDataSet'); 

set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')

figname = [dir_name, '/LOGsignaling2'];  
print(figname,'-depsc', '-r600'); 
% print(figname,'-dpng'); 
end

