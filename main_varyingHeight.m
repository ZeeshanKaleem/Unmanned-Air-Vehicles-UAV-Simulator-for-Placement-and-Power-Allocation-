%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Main Loop Runner in parallel:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear;
close all;
%clc;
NUAVs=8;
n=1;
ploting =0; % 0 to not plot, and 1 to plot
Q_learn =1; % Make it 1 for deciding actions based on Q-learning, 0 to allocat waterfilling based power.
scenario = 'k-means';%['k-means', 'fixed', 'random']; % Select a scenario of your choice
power_allocation_algo = 'PSO-PA'; % PSO-PA'Q_learning', 'Waterfilling','Max_power', 'Exhaustive Search','EPA','PSO' 
%iterations =100:200:5000;
iterations =50000;
rewards = 'distance'; %'double_rate', 'distance';
h = 1:40:700;
cresult = cell(1,n);
Result =cell(length(h),1);
ResultReward =zeros(1,length(iterations));
height_vs_rate = zeros(length(h),2);
f=2;
aplha=0.8;
for h1= 1:length(h)
  
for loops = 1:n
    %loops
saveNum = loops;
actions = 31; % We are considering 31 actions, that considers different power levels in each state. In this simulations we are focusing on power allocation only. Height and other dimenstions will remain static during the simulations.
permuted=1;
permutationsMat = randperm(NUAVs,NUAVs);


FBSSet_in = cell(1,1);


      
    [Q, C_Result] = PA_RL_permutatedUAVs(NUAVs,permutationsMat, 10,saveNum, h(h1),ploting,Q_learn,scenario,power_allocation_algo,rewards,h,h1,iterations);
    %[Q, C_Result] = PA_RL_permutatedUAVs(NUAVs,permutationsMat, 10,saveNum,100,ploting,Q_learn,scenario,power_allocation_algo,rewards,h,500,iterations,power(p)); 
    Q_out=[];
    Q_out= Q;
    result = [];
    result= C_Result;
    cresult{loops}= result;
end

%% Cumulative result struct


Result{h1} = cresult;

end
f_result = Result;
save(sprintf('results/fresult_%d_%d_%d_%s_%s_%s_%d.mat',NUAVs, length(h),n,scenario, power_allocation_algo, datetime('today'),h(h1)),'f_result')

