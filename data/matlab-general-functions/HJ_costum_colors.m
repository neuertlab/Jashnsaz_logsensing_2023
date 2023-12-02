function cmap = HJ_costum_colors(nn,dc)
%% dark

% dark 1 90min
col = [0 0 0]; 
for i=1:nn
    cmap{1}(i,:)=col+(1-col)*dc*(i-1);
end

% dark 1 35min
col = hex2rgb('1e2022'); 
col = [.2 .2 .2]; 
for i=1:nn
    cmap{2}(i,:)=col+(1-col)*dc*(i-1);
end

% dark 2 30min
col = hex2rgb('52616a'); 
col = [.4 .4 .4]; 
for i=1:nn
    cmap{3}(i,:)=col+(1-col)*dc*(i-1);
end

% dark 3 25min
col = hex2rgb('c9d6de'); 
col = [.6 .6 .6]; 

for i=1:nn
    cmap{4}(i,:)=col+(1-col)*dc*(i-1);
end

% dark 4 20 min
col = hex2rgb('f0f5f9');
col = [.8 .8 .8]; 

for i=1:nn
    cmap{5}(i,:)=col+(1-col)*dc*(i-1);
end

% teal 20 min 
col = hex2rgb('7f9eb2'); 
for i=1:nn
    cmap{6}(i,:)=col+(1-col)*dc*(i-1);
end

% brown  ctrl
col = hex2rgb('a7a7a2'); 
for i=1:nn
    cmap{7}(i,:)=col+(1-col)*dc*(i-1);
end
%% purple magenta 

% magenta
col = hex2rgb('9055A2'); 
for i=1:nn 
    cmap{8}(i,:)=col+(1-col)*dc*(i-1);
end

% purple
col = hex2rgb('791E94'); 
for i=1:nn 
    cmap{9}(i,:)=col+(1-col)*dc*(i-1);
end

%% blue cyan

% blue
col = hex2rgb('1500FE'); %[0 0 0.9]; 
for i=1:nn
    cmap{12}(i,:)=col+(1-col)*dc*(i-1);
end

% cyan
col = hex2rgb('2AC6FF');  % 84B1ED
for i=1:nn 
    cmap{11}(i,:)=col+(1-col)*dc*(i-1);
end

%% red orange yellow
% orange
col = hex2rgb('fc913a'); 
for i=1:nn 
    cmap{10}(i,:)=col+(1-col)*dc*(i-1);
end

% magenta
col = [1  0  1]; 
for i=1:nn
    cmap{13}(i,:)=col+(1-col)*dc*(i-1);
end

% copper
copp = copper(11); % copper, 
col = copp(4,:); % purple 791E94, 
for i=1:nn
    cmap{14}(i,:)=col+(1-col)*dc*(i-1);
end

% red
col = hex2rgb('c03546'); 
col = [1 0 0]; 
for i=1:nn
    cmap{15}(i,:)=col+(1-col)*dc*(i-1);
end

% yellow
col = hex2rgb('FFE955'); %EFDC05
for i=1:nn
    cmap{16}(i,:)=col+(1-col)*dc*(i-1);
end
end