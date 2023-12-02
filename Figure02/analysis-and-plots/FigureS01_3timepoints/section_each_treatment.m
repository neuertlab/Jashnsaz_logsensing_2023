function [timepoints, OD600] = section_each_treatment(t,OD)
% time  points
idd = diff([true,isnan(t),true]);
idb = find(idd<0);
ide = find(idd>0)-1;
timepoints = arrayfun(@(b,e)t(b:e),idb,ide,'uni',0); 

% OD600
idd = diff([true,isnan(OD),true]);
idb = find(idd<0);
ide = find(idd>0)-1;
OD600 = arrayfun(@(b,e)OD(b:e),idb,ide,'uni',0); 
end