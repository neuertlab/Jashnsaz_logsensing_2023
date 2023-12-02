function [Hog1SignalingData, dp] = get_BiolReps_stats(exp, dd, Hog1SignalingData)
    try
        n_BRs = size(Hog1SignalingData(exp).Sel,2); 
        for br=1:n_BRs
            RL(br) = Hog1SignalingData(exp).Sel(br).ttrans1(end);  % combine for loops
        end 
        
        Hog1m=NaN(max(RL),n_BRs); 
        Volm=NaN(max(RL),n_BRs);

        for br = 1:n_BRs
            Hog1m(1:RL(br),br) = nanmedian(Hog1SignalingData(exp).Sel(br).InIcSelNintpSBack2N(:,1:RL(br)),1); %nanmedian
            Volm(1:RL(br),br) = nanmedian(Hog1SignalingData(exp).Sel(br).areaN(:,1:RL(br)),1);
            tt(1:RL(br),br) = Hog1SignalingData(exp).Sel(br).ttrans1(1:RL(br))-Hog1SignalingData(exp).Sel(br).t0;
        end
        
        H = NaN(max(RL(:))+50,n_BRs); V = H; t = H; t1 = t; H1 = H; V1 = V; sRR = n_BRs; %dp = V;
        
        for i = 1:n_BRs            
            t(1:RL(i),i) = tt(1:RL(i),i)';
            H(1:RL(i),i) = Hog1m(1:RL(i),i)';
            V(1:RL(i),i) = Volm(1:RL(i),i)';
        end
        
        % Correct for time offset between experiemtns
        
%         for i = 1:sRR;       
            ts1 = t(1,:);
            [ts2,its2] = max(ts1);
            ts3 = ts1-ts2;
%         end
        
        for i = 1:sRR;                               
            t1(:,i) = circshift(t(:,i),ts3(i));
            H1(:,i) = circshift(H(:,i),ts3(i));
            V1(:,i) = circshift(V(:,i),ts3(i));
        end
%         size(V1')
%         size(H1')
        % Determine mean from replica experiments
        Hog1SignalingData(exp).Hog1m = nanmean(H1',1); 
        Hog1SignalingData(exp).Hog1s = nanstd(H1',0,1);
        Hog1SignalingData(exp).Volm = nanmean(V1',1);
        Hog1SignalingData(exp).Vols = nanstd(V1',0,1);
        Hog1SignalingData(exp).tt = nanmean(t1',1)./6;
                
%         maxID = min([size(Hog1SignalingData(exp).Hog1m,2), size(Hog1SignalingData(exp).Hog1s,2), size(Hog1SignalingData(exp).Volm,2), size(Hog1SignalingData(exp).Vols,2), size(Hog1SignalingData(exp).tt,2)]);  
%         Hog1SignalingData(exp).tt = Hog1SignalingData(exp).tt(1:maxID);
%         Hog1SignalingData(exp).Hog1m = Hog1SignalingData(exp).Hog1m(1:maxID);
%         Hog1SignalingData(exp).Hog1s = Hog1SignalingData(exp).Hog1s(1:maxID);        
%         Hog1SignalingData(exp).Volm = Hog1SignalingData(exp).Volm(1:maxID);
%         Hog1SignalingData(exp).Vols = Hog1SignalingData(exp).Vols(1:maxID); 

        for i=1:size(Hog1SignalingData(exp).tt,2)
            if i> size(Hog1SignalingData(exp).tt,2)/2 & Hog1SignalingData(exp).tt(i)<=0
                Hog1SignalingData(exp).tt(i) = NaN;
            end
        end
    end  

        AA = floor(Hog1SignalingData(exp).tt) == Hog1SignalingData(exp).tt;
        dpA = find(AA == 1,1);
        dpE = find(AA == 1,1,'last');
        dp = [dpA:dd:dpE]; 
        
        Hog1SignalingData(exp).dp = dp;

end