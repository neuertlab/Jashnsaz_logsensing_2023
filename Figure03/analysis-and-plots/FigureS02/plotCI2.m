function plotCI2(xvar, var, col)
    ploting = 1 ; % col=ones(1,3); 
%     fill_between_lines1 = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C,'facealpha',.1,'edgecolor','none');
    fill_between_lines2 = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C,'facealpha',.3,'edgecolor',col);

%     fill_between_lines1(xvar,prctile(var, 30), prctile(var, 70), col); hold on; 
    fill_between_lines2(xvar,prctile(var, 20), prctile(var, 60), col); hold on; 

end 
