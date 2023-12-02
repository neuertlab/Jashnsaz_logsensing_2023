function [cmap] = colormapHJ(N,hex_id)
    cmap=NaN(N,3); dd=(1/(1.3*N)); 
    col = hex2rgb(hex_id);
    for i=1:N
        cmap(i,:)=col+(1-col)*dd*(i-1);
    end
end