%% Plots
% For 4 UAVs case and n loops, height lenght

clear;
h = 100:400:800;
NUAVs =4; h_lenght = length(h); n=2;
s = sprintf('fresult_%d_%d_%d.mat',NUAVs,h_lenght,n);
filename = strcat(s);
load(filename);
for h1= 1:length(h) 

% Sum rate 
res =f_result{h1,:};
r=0;
sum =0;

min_cuav = 0;
sum_minCFUE=0;
for i=1:n
r = res{i}.sum_CFUE;
sum = sum +r;

min_cuav = res{i}.min_CFUE;
sum_minCFUE = sum_minCFUE +min_cuav;
end


mean_sum_rate = sum/n;
mean_cuav_rate = sum_minCFUE/n;

height_vs_rate(h1,:) = [h(h1) mean_sum_rate];
height_vs_uavrate(h1,:) = [h(h1) mean_cuav_rate];
end
%% plot height vs mean sum rate
figure
plot(height_vs_rate(:,1), height_vs_rate(:,2), 'r*');
%% plot height vs UAV UE rate
figure
plot(height_vs_uavrate(:,1), height_vs_uavrate(:,2), 'k-');
%end
