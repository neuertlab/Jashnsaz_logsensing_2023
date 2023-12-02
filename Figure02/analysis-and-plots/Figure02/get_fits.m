function [t0,f,CI_f,mean_T,neg_CI,pos_CI] = get_fits(t,OD600,options)
% define fit options and fit type. 
if options.drop=="yes"
    t(1:options.howmany)=[]; 
    OD600(1:options.howmany)=[]; 
end
t0=t(1); 
t=t-t0;
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...                                             %[doubling_time] > [0]
               'Upper',[inf],...                                           %[doubling_time] < [inf]
               'StartPoint',[3000+(2000*rand-1000)]);                      % guess for [doubling_time], uniform distribution

ft = fittype('OD0*2^(x/b)','problem','OD0','options',fo);                   % define fit type: 2^(t/T) growth rate function OD0 is OD600 value at t=0). 

% fit OD600 measurment (t, OD600) to an exponential function. 
[f] = fit(t,OD600,ft,'problem',OD600(1));
CI_f = confint(f,0.95);

mean_T=[f.b]; 
neg_CI=CI_f(1)-f.b; 
pos_CI=CI_f(2)-f.b;
end
% figure;histogram([3000+(2000*rand(1000,1)-1000)])