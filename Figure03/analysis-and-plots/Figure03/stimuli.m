ExperimentsID = {[1:6], [7:12], [13:18], [19:24], [25:30], [31:36], [37:42]}; 
Experiments = {'steps $ (t^{0}) $', 'root $ (\sqrt{t}) $', 'linears $ (t^{1}) $', 'quadratics $ (t^{2}) $','cupics $ (t^{3}) $', 'quintics $ (t^{5}) $', 'heptics $ (t^{7}) $'}; 
dir_name = 'figures'; mkdir(dir_name);
close all

% set total time and NaCl on time
TTT = 50; TT=25; 

% pick functions from roots and polynomials
roots = [2]; 
polynomials = [1 2 3 5 7]; 

% Final NaCl molarities (M)
NaCls = [0 .1 .2 .4 .6 .8]; 

n_polynomials = length(polynomials); 
n_roots = length(roots);
n_NaCls = length(NaCls); 
all_funs = n_polynomials+n_roots; 
tot_inputs = all_funs*n_NaCls; 

fxngrps = 1; 
% step
sympref('HeavisideAtOrigin',1);
step_fun = @(t, maxNaCl) maxNaCl * heaviside(t - .001);

fxngrps = fxngrps + 1; 
% root functions
for i=1:n_roots 
    root_fun{i} = @(t, maxNaCl) min([(maxNaCl.*t.^(1/roots(i))/(60*TT)^(1/roots(i))),maxNaCl]); 
    fxngrps = fxngrps + 1; 
end

% polynomial functions
for i=polynomials
    poly_fun{i} = @(t, maxNaCl) min([(maxNaCl.*t.^i/(60*TT)^i),maxNaCl]);
    fxngrps = fxngrps + 1; 
end
% power law function inputs
pow_inputs = @(t, maxNaCl) min([(-1 + ((1+maxNaCl)^(1/(60*TT)))^t),maxNaCl]);

fun_n = 1;
for steps=NaCls
    stimuli_funcs{fun_n} = @(t) step_fun(t, steps); 
    fun_n = fun_n +1; 
end

for rootfuns = 1:n_roots
    for ramps=NaCls
        stimuli_funcs{fun_n} = @(t) root_fun{rootfuns}(t, ramps); 
        fun_n = fun_n +1; 
    end
end

for polys = polynomials
    for ramps=NaCls
        stimuli_funcs{fun_n} = @(t) poly_fun{polys}(t, ramps); 
        fun_n = fun_n +1; 
    end
end

for steps=NaCls
    stimuli_funcs{fun_n} = @(t) pow_inputs(t, steps); 
    fun_n = fun_n +1; 
end
fun_n
for i=1:fun_n-1
    data{i}.tt=[-5:.1:57];
    data{i}.stimuli = stimuli_funcs{i};
    
    j=1;
    stimuli_profile=zeros(1,length(data{i}.tt));
    for t = data{i}.tt
        if t<0
            stimuli_profile(j) = 0; 
        else
            stimuli_profile(j) = data{i}.stimuli(t*60);
        end
        j=j+1;
    end 
    data{i}.stimuli_profile=stimuli_profile;  
end

% %% plot stimuli
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
%         plot(data{exp}.tt,data{exp}.stimuli_profile,'LineWidth', lw, 'color',cmap(count,:)); 
%         box on; xlim([-.5 50.5]); ylim([-.01 .81]); xticks([0:10:50]); yticks([0:.2:.8]); 
%         count = count + 1; 
%     end
% end
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',8, 'defaultTextFontSize',8, 'FontName', 'Helvetica')
% figname = [dir_name, '/NaCl'];  
% print(figname,'-depsc', '-r600');
%%
close; figure(2); set(gcf, 'Units', 'centimeters', 'Position', [0 0 21.5 6], 'PaperUnits', 'centimeters', 'PaperSize', [21.5 8]); 
dx = 0.05; dy = 0.08; 

cmap = hsv(7);

cmap = [0 0 0; cmap]; 
count = 2; 

subplotHJ(1,4,1,dy,dx); hold on; grid on;
set(groot,'defaultAxesTickLabelInterpreter','tex');

plot(data{1}.tt,zeros(size(data{1}.tt)),'LineWidth', lw, 'color',cmap(1,:));

for i=1:size(ExperimentsID,2)    
        
    for exp = ExperimentsID{i}(end-1)
        plot(data{exp}.tt,data{exp}.stimuli_profile,'LineWidth', lw, 'color',cmap(count,:)); 
        box on; xlim([-.5 50.5]); ylim([-.01 .81]); xticks([0:10:50]); yticks([0:.2:.8]); 
        count = count + 1; 
    end
end
box on; xlim([-3.5 50.5]); ylim([-.01 .61]); xticks([0:10:50]); yticks([0:.2:.6]); 
