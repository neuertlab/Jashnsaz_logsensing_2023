function [data] = stimuli_functions()
% set total time and NaCl on time
TTT = 50; TT=25; 

% pick functions from roots and polynomials
roots = [2]; 
polynomials = [1 2 3 5 7]; 

% Final NaCl molarities (M)
NaCls = [0 .1 .2 .3 .4 .6 .8]; 

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

for i=1:fun_n-1
    data{i}.tt=[-30:.1:57];
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

end

