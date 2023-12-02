function plotCI(xvar, var, col)
    ploting = 1 ; % col=ones(1,3); 
    fill_between_lines1 = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C,'facealpha',.1);
    fill_between_lines2 = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C,'facealpha',.4);

    fill_between_lines1(xvar,prctile(var, 40), prctile(var, 60), col); hold on; 
    fill_between_lines2(xvar,prctile(var, 45), prctile(var, 55), col); hold on; 

end 
