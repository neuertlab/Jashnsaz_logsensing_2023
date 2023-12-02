function [error, best_parameters, Vol_fit] = polynomial_function_fit(fxncall,data,par)
    
    %% parameters
        % y0, delta, [0,1]
        % k, polynomial order, [.5,1,2,3, 7]
        % T0, duration, 25min
        % T1, delay, 1min 
        % T2, expansion, 1min
%         delta = [.5 .1 1 .1 .1];
        delta = [.5 0 0 0 0]; 
        if fxncall==2; delta(2)=3; end
        if fxncall==3; delta([2 5])=[3 2]; end

    %% optimize
    n_params = numel(delta); 
    initial_parameters.k = [.5 par(1) par(2) 0 0];
    % set optimization
    parameters_limits = optimvar('k',n_params,"LowerBound",initial_parameters.k-delta...
        ,"UpperBound",initial_parameters.k+delta);

    obj = 0;
    model = fcn2optimexpr(@polynomial_fxn,data.time,parameters_limits);
    obj = obj+sum(sum((model - data.vol).^2));

    problem = optimproblem("Objective",obj);

    % best_parameters_chain  = []; 
    initial_parameters.k = initial_parameters.k; % +1e-2*ones(1,n_params);
    if fxncall==2; initial_parameters.k(2)=10*rand; end
    if fxncall==3; initial_parameters.k([2 5])=[10*rand 1]; end
%     initial_parameters.k = [.5 10*rand 25+randn randn randn];

    [best_parameters,error] = solve(problem,initial_parameters);
    % best_parameters_chain=[best_parameters_chain; best_parameters.k']; 
    Vol_fit = polynomial_fxn(data.time,best_parameters.k); 
    best_parameters = best_parameters.k; 
end


