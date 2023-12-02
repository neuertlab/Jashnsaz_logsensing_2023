function sol = polynomial_fxn(t,param)

%     if jcall == 2 % step
%         polynomial_function  = @(t,y0,n,T0,T1,T2)      [0].*(t<=T1)...
%                                                     + [y0].*(t>T1);
%     else
    polynomial_function = @(t,y0,n,T0,T1,T2)       [0].*(t<=T1)...
                                                + [y0*((t-T1)./(T0-T1+T2)).^n].*((t>T1)&(t<=T0+T1+T2))...
                                                + [y0].*(t>T0+T1+T2);
%     end
    sol = polynomial_function(t,param(1),param(2),param(3),param(4),param(5));
end

    % true_params = [.5 25 0 0.5 2];
    % time = (-10:.1:50)';
    % 
    % % figure(); hold on
    % plot(time,polynomial_fxn(time,true_params))