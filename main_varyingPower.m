%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Main Loop Runner in parallel:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear;
close all;
%clc;
NUAVs=4;
n=1;
ploting =0; % 0 to not plot, and 1 to plot
Q_learn =1; % Make it 1 for deciding actions based on Q-learning, 0 to allocat waterfilling based power.
scenario = 'k-means';%['k-means', 'fixed', 'random']; % Select a scenario of your choice
power_allocation_algo = 'EPA'; % PSO-PA'Q_learning', 'Waterfilling','Max_power', 'Exhaustive Search','EPA','PSO' 
%iterations =100:200:5000;
iterations =50000;
rewards = 'distance'; %'double_rate', 'distance';
h = 1:40:700;
cresult = cell(1,n);
%Result =cell(length(h),1);
ResultReward =zeros(1,length(iterations));
height_vs_rate = zeros(length(h),2);
f=2;
aplha=0.8;
power=20:20:300;
Result =cell(length(power),1);
%for iter= 1:length(iterations) % Activate when have to plot the reward
%curve
for p= 1:length(power)
  
for loops = 1:n
    %loops
saveNum = loops;
%states = 28; % We are only using four fixed states
actions = 31; % We are considering 31 actions, that considers different power levels in each state. In this simulations we are focusing on power allocation only. Height and other dimenstions will remain static during the simulations.
permuted=1;
%QTable= zeros (NUAVs,actions);

%permutationsMat = zeros(100,16);
%permutationsMat = zeros(1,NUAVs);

permutationsMat = randperm(NUAVs,NUAVs);


FBSSet_in = cell(1,1);

%for i=1:16
      
    %[Q, C_Result] = PA_RL_permutatedUAVs(NUAVs,permutationsMat, 10,saveNum, h(h1),ploting,Q_learn,scenario,power_allocation_algo,rewards,h,h1,iterations);
    [Q, C_Result] = PA_RL_permutatedUAVs(NUAVs,permutationsMat, 10,saveNum,100,ploting,Q_learn,scenario,power_allocation_algo,rewards,500,500,iterations,power(p)); 
    Q_out=[];
    Q_out= Q;
    result = [];
    result= C_Result;
    cresult{loops}= result;

end

%% Cumulative result struct
%c_answer.result = cresult;

Result{p} = cresult;

end
%ResultReward(:,iter) =result.AvgR; 
%end
f_result = Result;
save(sprintf('results/fresult_%d_%d_%d_%s_%s_%s_%d.mat',NUAVs, length(h),n,scenario, power_allocation_algo, datetime('today'),power(p)),'f_result')
