function [yyy,fit_fun] = expfit(xx0,yy)

% xx  = [xx0 100 120 200]; yy=[yy yy(end)-(1e-3*rand) yy(end)-(1e-2*rand) yy(end)-(1e-2*rand)]; 
xx  = [xx0 linspace(xx0(end),200,10)]; yy=[yy linspace(yy(end), yy(end)/2,10)]; 

ff = fit(xx',yy','exp2');                                               % fit a double exponential function as background
fit_fun = @(t) ff.a*exp(ff.b*t) + ff.c*exp(ff.d*t); 
yyy = fit_fun(xx0); 

end

