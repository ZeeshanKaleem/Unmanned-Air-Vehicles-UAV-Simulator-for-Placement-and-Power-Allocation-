%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of Power Allocation in femtocell network using 
%   Reinforcement Learning with random adding of femtocells to the network
%   
%
function [Q, QFinal] = PA_RL_permutatedUAVs(fbsCount,femtocellPermutation, NumRealization, saveNum, h1, ploting, Q_learn,scenario,power_allocation_algo,rewards,h2,h3,iter,pmax)

%% Initialization
% clear all;
clc;
format short
format compact
tic;
%% Parameters
Pmin = -20;                                                                                                                                                                                                                                                                                                                                                                           %dBm
Pmax = 25; %dBm
%Pmax = 10*log10(pmax); %dBm
Npower = 31;
nMBS =1;
nMUE = 1;
dth = 25;
Kp = 100; % penalty constant for MUE capacity threshold
Gmue = 1.37; % bps/Hz
StepSize = 1.5; % dB
K = 1000;
PBS = 50 ; %dBm
noise= -174; % dBm
rad_mbs =500;
num_points=100;
range =500;
dth = 120; %meter % threshold distance to check how far the FBS is from neighbor FUE
power = zeros(fbsCount,fbsCount);
powervect =[];
sumratef=[];
powervectf = [];
small_scale =0;

%% UAV specific parameters
h= repelem(h1,fbsCount); % Fixed height 
varying_threshold=0;
rad_Ubs = 100;
switch scenario % This is valid for 4 UAVs only
    case 'fixed'
          
        if fbsCount==4
          xUAV= [range/2 -range/2, -range/2, range/2];   % We can Generate values from the uniform distribution on the interval (a, b) as: r = a + (b-a).*rand(5,1);
%        
          yUAV =[range/2 range/2, -range/2, -range/2];
          
        elseif fbsCount==8
            range = range/2;
                xUAV= [range/2; (range+range+range)/2; -range/2; (-range-range-range)/2; ...
                    -range/2; (-range-range-range)/2; range/2; (range+range+range)/2;];
    
                 yUAV =[(range+range)/2;(range+range)/2;  (range+range)/2; (range+range)/2; ...
                 (-range-range)/2;(-range-range)/2; (-range-range)/2; (-range-range)/2;];
              range =2*range; % this is for value 
        elseif fbsCount==16
            
            range = range/4;
                r=range;% 25
                r1=r+range; % 50
                r2 = r1+range; % 75
                r3 = r2+range; % 100


                xUAV= [r/2; (r+r1)/2;(r1+r2)/2; (r2+r3)/2;...
                    -r/2; (-r1-r)/2;(-r2-r1)/2; (-r3-r2)/2;...
                    -r/2; (-r1-r)/2;(-r2-r1)/2; (-r3-r2)/2;...
                    r/2; (r+r1)/2;(r1+r2)/2; (r2+r3)/2];

                yUAV =[4*r/2;4*r/2; 4*r/2; 4*r/2; ...
                    4*r/2; 4*r/2;4*r/2;4*r/2;...
                    -4*r/2;-4*r/2;-4*r/2;-4*r/2;...
                    -4*r/2;-4*r/2;-4*r/2;-4*r/2];
                 range = 4*range;
            
            
        end
    case 'k-means'

            [xUAV,yUAV] = optimal_x_y_points(num_points,range,fbsCount, ploting); % Optimal points calculation using k-means algorithm
            dis = [xUAV;yUAV];
    
    case 'random'
        
          if fbsCount==4
            xUAV= [randi(range,1,1), randi([-range 0],1,1), randi([-range 0],1,1), randi(range,1,1)];   % We can Generate values from the uniform distribution on the interval (a, b) as: r = a + (b-a).*rand(5,1);
       
            yUAV =[randi(range,1,1),  randi(range,1,1), randi([-range 0],1,1), randi([-range 0],1,1)];
          elseif fbsCount==8
               range = range/2;
                 xUAV= [randi([1 range],1,1); randi([range range+range],1,1); randi([-range 1],1,1); randi([-range-range -range],1,1); ...
                    randi([-range 1],1,1); randi([-range-range -range],1,1); randi([1 range],1,1);randi([range range+range],1,1)];
    
                 yUAV =[randi([1 range+range],1,1);randi([1 range+range],1,1);  randi([1 range+range],1,1); randi([1 range+range],1,1); ...
                 randi([-range-range 1],1,1);randi([-range-range 1],1,1); randi([-range-range 1],1,1); randi([-range-range 1],1,1)];
              range =2*range; % this is for value 
          elseif  fbsCount==16
                range = range/4;
                r=range;% 25
                r1=r+range; % 50
                r2 = r1+range; % 75
                r3 = r2+range; % 100


                xUAV= [randi([1 r],1,1); randi([r r1],1,1);randi([r1 r2],1,1); randi([r2 r3],1,1);...
                    randi([-r 1],1,1); randi([-r1 -r],1,1);randi([-r2 -r1],1,1); randi([-r3 -r2],1,1);...
                    randi([-r 1],1,1); randi([-r1 -r],1,1);randi([-r2 -r1],1,1); randi([-r3 -r2],1,1);...
                    randi([1 r],1,1); randi([r r1],1,1);randi([r1 r2],1,1); randi([r2 r3],1,1)];

                yUAV =[randi([1 4*r],1,1);randi([1 4*r],1,1);  randi([1 4*r],1,1); randi([1 4*r],1,1); ...
                    randi([1 4*r],1,1);randi([1 4*r],1,1);  randi([1 4*r],1,1); randi([1 4*r],1,1);...
                    randi([-4*r 1],1,1);randi([-4*r 1],1,1);  randi([-4*r 1],1,1); randi([-4*r 1],1,1);...
                    randi([-4*r 1],1,1);randi([-4*r 1],1,1);  randi([-4*r 1],1,1); randi([-4*r 1],1,1)];
                 range = 4*range;
            end
              
        end
%% Minimum Rate Requirements for N MUE users
N = 3;
q_mue = 1.00; q_fue=0.4;
step =2;

%% Q-Learning variables
% Actions
switch power_allocation_algo
    case 'Q_learning'
    actions= Pmin:step:Pmax;
    case 'Exhaustive Search'
   
    actions= Pmin:step:Pmax;
    otherwise
     actions= Pmin:step:Pmax;
end
Powers =zeros(length(actions),fbsCount);
        for i=1:length(actions)
                 Powers(i,:) = repelem(actions(i),fbsCount);
        end
 
%%
% to implement different combinations in exhaustive search
switch power_allocation_algo
    case 'Exhaustive Search'
   if fbsCount==4
   
    Combin = allcomb(actions,actions,actions,actions);
elseif fbsCount==8
    Combin = allcomb(actions,actions,actions,actions,actions,actions,actions,actions); 
else
     Combin = allcomb(actions,actions,actions,actions,actions,actions,actions,actions,actions,actions,actions,actions,actions,actions,actions,actions); 
   end
end
%Combin = nchoosek(actions,fbsCount);
%%
% States
%states = allcomb(0:3 , 0:3); % states = (dMUE , dBS)
if fbsCount==4
    states = allcomb(0:1 , 0:1);
elseif fbsCount==8
    states = allcomb(0:1 , 0:1, 0:1); % states for fixed x,y position and varying Height.
else
    states = allcomb(0:1 , 0:1, 0:1, 0:1);
end


% Q-Table
% Q = zeros(size(states,1) , size(actions , 2));
Q1 = ones(size(states,1) , size(actions , 2)) * inf;
QTable= zeros (size(states,1),size(actions , 2));


alpha = 0.9; gamma = 0.5; epsilon = 0.1; 
switch power_allocation_algo
    case 'Exhaustive Search'
    Iterations = size((Combin),1);
    case 'Q_learning'
    Iterations = iter;
    otherwise
    Iterations = 100;    
end
%% Generate the UEs
 %mue(1) = UE(204, 207);
% mue(1) = UE(150, 150);
% mue(1) = UE(-200, 0);
% selectedMUE = mue(mueNumber);
MBS = BaseStation(0 , 0 , 50);

x=0;
for i=1:nMBS
    xMUE(1+x:nMUE+x)=(MBS.X-rad_mbs) + ((MBS.X+rad_mbs)-(MBS.X-rad_mbs))*rand(nMUE,1);
    yMUE(1+x:nMUE+x)=(MBS.Y-rad_mbs) + ((MBS.Y+rad_mbs)-(MBS.Y-rad_mbs))*rand(nMUE,1); 
    x=x+nMUE;
end

for i=1:nMUE
    mue(i) = UE(xMUE(i), yMUE(i));
end
%%
%Generate fbsCount=16 FBSs
FBS_Max = cell(1,fbsCount);

for i=1:fbsCount
%     if i<= fbsCount
        FBS_Max{i} = FemtoStation_3S(xUAV(i),yUAV(i), h(i), MBS, mue, 10,rad_Ubs);
end
% for i=1:3
%     if i<= fbsCount
%         FBS_Max{i} = FemtoStation_3S(180+(i-1)*35,150, h, MBS, mue, 10,rad_Ubs);
%     end
% end
% 
% for i=1:3
%     if i+3<= fbsCount
%         FBS_Max{i+3} = FemtoStation_3S(165+(i-1)*30,180, h, MBS, mue, 10,rad_Ubs);
%     end
% end
% 
% for i=1:4
%     if i+6<= fbsCount
%         FBS_Max{i+6} = FemtoStation_3S(150+(i-1)*35,200,h, MBS, mue,10,rad_Ubs);
%     end
% end
% 
% for i=1:3
%     if i+10<= fbsCount
%         FBS_Max{i+10} = FemtoStation_3S(160+(i-1)*35,240, h, MBS, mue,10,rad_Ubs);
%     end
% end
% 
% for i=1:3
%     if i+13<= fbsCount
%         FBS_Max{i+13} = FemtoStation_3S(150+(i-1)*35,280, h, MBS, mue, 10,rad_Ubs);
%     end
% end


%%
% 
FBS = cell(1,fbsCount);

temp = [];
temp2 =[];
temp3= [];
temp4= [];
for i=1:fbsCount
    FBS{i} = FBS_Max{femtocellPermutation(i)};
    x = FBS{i}.X;
    xUAV= [temp x];
    temp =xUAV;
    
    y = FBS{i}.Y;
    yUAV= [temp2 y];
    temp2 =yUAV;
    
    fuex= FBS{i}.FUEX;
    xUvUE =[temp3 fuex];
    temp3 =xUvUE;
    
    fuey= FBS{i}.FUEY;
    yUvUE =[temp4 fuey];
    temp4 =yUvUE;
end

%% Plotting System Model
if ploting
figure;
% plot(MBS.X,MBS.Y,'>k','Linewidth',2)
% hold on;
plot(xUAV,yUAV,'>r','Linewidth',2)
hold on;
% plot(xMUE,yMUE,'sr','Linewidth',1)
% hold on;
plot(xUvUE,yUvUE,'og','Linewidth',1)


% circle(MBS.X,MBS.Y,rad_mbs, 'r');
circle(FBS{1}.X,FBS{1}.Y,rad_Ubs, 'g');
circle(FBS{2}.X,FBS{2}.Y,rad_Ubs, 'k');
circle(FBS{3}.X,FBS{3}.Y,rad_Ubs, 'b');
circle(FBS{4}.X,FBS{4}.Y,rad_Ubs, 'r');
circle(FBS{5}.X,FBS{5}.Y,rad_Ubs, 'g');
circle(FBS{6}.X,FBS{6}.Y,rad_Ubs, 'k');
circle(FBS{7}.X,FBS{7}.Y,rad_Ubs, 'b');
circle(FBS{8}.X,FBS{8}.Y,rad_Ubs, 'r');
% circle(selectedMUE.X,selectedMUE.Y,dM2, 'r');
% circle(selectedMUE.X,selectedMUE.Y,dM3, 'r');
end
    %% Initialization and find MUE Capacity
    % permutedPowers = npermutek(actions,3);
    switch power_allocation_algo
        case'Exhaustive Search'
        
        %textprogressbar((episode/Iterations)*100);
%          step_size = 46/size(FBS,2);
%         permutedPowers = -20:step_size:25;
%         
%         
%         for i=1:length(permutedPowers)
%             for j=1:length(permutedPowers)
%                 power(i,j)=permutedPowers(j);
%             end
%         end
%         horzon = horzcat(permutedPowers',power);
%         ee=[];
%         for jj=1:size(horzon,1)
%             cc= horzon(jj,:);
%             dd= unique(cc,'stable');
%             ee= [ee;dd];
%         end
          %count=count+1;
          %c= randi(size(FBS,2),1);
          c= randi(size(Combin,1));
            % Action selection with epsilon=0.1
            %if rand<epsilon   % For implementation of epsilon greedy algorithm (It will explore with the probability of epsilon) : Exploration here means that randomly search for other power options
                for j=1:size(FBS,2)
                    fbs = FBS{j};
                    fbs = fbs.setPower(Combin(c,j));
                    FBS{j} = fbs;
                end
            
        
       % a=repelem(actions,size(FBS,2));
       % b=reshape(a,[size(FBS,2),size(actions,2)]);
        %permutedPowers = randperm(size(actions,2),size(FBS,2));
        %permutedPowers = randi(size(actions,2)*size(FBS,2),1, size(FBS,2));
    
        for j=1:size(FBS,2)
            fbs = FBS{j};
            %fbs = fbs.setPower(b(permutedPowers(j)));
            fbs = fbs.getDistanceStatus(range, fbsCount,scenario);
            d_min_neigh = nearest_FUE(fbs,FBS,fbsCount,j);
            fbs = fbs.min_neighbordistance(d_min_neigh);
            FBS{j} = fbs;
        end
     case'Q_learning'
        permutedPowers = randperm(size(actions,2),size(FBS,2));
    
    % y=randperm(size(permutedPowers,1));
        for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setPower(actions(permutedPowers(j)));
            fbs = fbs.getDistanceStatus(range, fbsCount,scenario);
            d_min_neigh = nearest_FUE(fbs,FBS,fbsCount,j);
            fbs = fbs.min_neighbordistance(d_min_neigh);
            FBS{j} = fbs;
        end
        
        case'EPA'
            sel = randperm(length(actions),1);
            sel_p = actions(sel);
        Powersp = repelem(sel_p,fbsCount);
    
    % y=randperm(size(permutedPowers,1));
        for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setPower(Powersp(j));
            fbs = fbs.getDistanceStatus(range, fbsCount,scenario);
            d_min_neigh = nearest_FUE(fbs,FBS,fbsCount,j);
            fbs = fbs.min_neighbordistance(d_min_neigh);
            FBS{j} = fbs;
        end
        
        case 'PSO'
            
%             for j=1:size(FBS,2)
%                     fbs = FBS{j};
%                     fbs = fbs.setPower(Pmax);
%                     FBS{j} = fbs;
%             end
%             
            
            sel = randperm(length(actions),1);
            sel_p = actions(sel);
            Powersp = repelem(sel_p,fbsCount);
             for j=1:size(FBS,2)
                    fbs = FBS{j};
                   fbs = fbs.setPower(Powersp(j));
                    FBS{j} = fbs;
             end
        
        
            
           case 'PSO-PA' 
              sel = randperm(length(actions),1); %intially assign equal power among the base stations
             sel_p = actions(sel);
            Powersp = repelem(sel_p,fbsCount);
             for j=1:size(FBS,2)
                    fbs = FBS{j};
                   fbs = fbs.setPower(Powersp(j));
                    FBS{j} = fbs;
             end
            
        otherwise
            
            permutedPowers = randperm(size(actions,2),size(FBS,2));
    
    % y=randperm(size(permutedPowers,1));
        for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setPower(actions(permutedPowers(j)));
            fbs = fbs.getDistanceStatus(range, fbsCount,scenario);
            d_min_neigh = nearest_FUE(fbs,FBS,fbsCount,j);
            fbs = fbs.min_neighbordistance(d_min_neigh);
            FBS{j} = fbs;
        end
        
                  
    end
%     selectedMUE.SINR = SINR_MUE(FBS, BS, selectedMUE, -120, 1000);
%     selectedMUE.C = log2(1+selectedMUE.SINR);

%     if selectedMUE.C < gamma_th
%         I = 1;
%     else
%         I = 0;
%     end
% 
%     for j=1:size(FBS,2)
%         fbs = FBS{j};
%         fbs.state(1,1) = I;
%         FBS{j} = fbs;
%     end
%% Calc channel coefficients
    fbsNum = size(FBS,2);
    if size(mue,2)<2
        G = zeros(fbsNum+1, fbsNum+1); % Matrix Containing small scale fading coefficients
        L = zeros(fbsNum+1, fbsNum+1); % Matrix Containing large scale fading coefficients
    else
        G = zeros(fbsNum+1, fbsNum+size(mue,2)); % Matrix Containing small scale fading coefficients
        L = zeros(fbsNum+1, fbsNum+size(mue,2)); % Matrix Containing large scale fading coefficients
    end
        
    
       switch power_allocation_algo
           case 'PSO'
             positions=PSO(FBS,MBS,mue,NumRealization,fbsCount,h,10,rad_Ubs);
                
              ii=1;
              jj=2;
              for a=1:size(FBS,2)
                  fbs = FBS{a};
                  fbs.X=positions(:,ii);
                  fbs.Y=positions(:,jj);
                    ii=ii+2;
                    jj=jj+2;
                  
                  
                fbs1 = FemtoStation_3S(fbs.X,fbs.Y, h, MBS, mue, 10,rad_Ubs);
                fbs.FUEX = fbs1.FUEX;
                fbs.FUEY = fbs1.FUEY;
                FBS{a}=fbs;

              end
              [G,L]= measure_channel_UAV_PSO1(FBS,MBS,mue,NumRealization, varying_threshold);
          otherwise
           [G, L] = measure_channel_UAV(FBS,MBS,mue,NumRealization, varying_threshold);
       end
       
       
    %% Main Loop
    fprintf('Loop for %d number of FBS :\t', fbsCount);
    textprogressbar(sprintf('calculating outputs:'));
    count = 1;
    MUE_C = zeros(1,Iterations);
    xx = zeros(1,Iterations);
    errorVector = zeros(1,Iterations);
    accumulated_reward=zeros(1,Iterations);
    temp_sumrate = zeros(1,Iterations);
    powervector = zeros(Iterations,fbsCount);
        % K1 is distance of selectedMUE from Agents
%     k1 = zeros(1,size(FBS,2));
    
%     Kp = 100;
%     for i=1:size(FBS,2)
%         k1(i) = (sqrt((FBS{i}.X-selectedMUE.X)^2+(FBS{i}.Y-selectedMUE.Y)^2))/dth;
%     end
    for episode = 1:Iterations
       switch power_allocation_algo
            case 'Q_learning'
        textprogressbar((episode/Iterations)*100);
        permutedPowers = randperm(size(actions,2),size(FBS,2));
        if (episode/Iterations)*100 < 80
            % Action selection with epsilon=0.1
            if rand<epsilon   % For implementation of epsilon greedy algorithm (It will explore with the probability of epsilon) : Exploration here means that randomly search for other power options
                for j=1:size(FBS,2)
                    fbs = FBS{j};
                    fbs = fbs.setPower(actions(permutedPowers(j)));
                    FBS{j} = fbs;
                end
            else
                for j=1:size(FBS,2)
                    fbs = FBS{j};
                    for kk = 1:size(states,1)
                        if states(kk,:) == fbs.state
                            break;
                        end
                    end
                    [M, index] = max(QTable(kk,:));
                    fbs = fbs.setPower(actions(index));
                    FBS{j} = fbs;
                end
            end
        else % It means if 80% of iterations are already being performed, then don't need to explore more.
            for j=1:size(FBS,2)
                fbs = FBS{j};
                for kk = 1:size(states,1)
                    if states(kk,:) == fbs.state
                        break;
                    end
                end
                [M, index] = max(QTable(kk,:));
                fbs = fbs.setPower(actions(index));
                FBS{j} = fbs;
            end
        end 

        % calc FUEs and MUEs capacity
        SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, noise,small_scale);
        C_FUE_Vec = log2(1+SINR_FUE_Vec);
        for i=1:size(mue,2)
            MUE = mue(i);
            MUE.SINR = SINR_MUE_4(G, L, FBS, MBS, MUE, noise);
            MUE = MUE.setCapacity(log2(1+MUE.SINR));
            mue(i)=MUE;
        end
        
%         MUE_C(1,episode) = selectedMUE.C;
        xx(1,episode) = episode;
%         R = K - (selectedMUE.SINR - sinr_th)^2;
%             deviation_FUE=0.0;
%             for i=1:size(FBS,2)
%                 deviation_FUE = deviation_FUE + (fbs.C_FUE-q_M)^2;
%             end
%         dum1 = 1.0;
%         for i=1:size(mue,2)
%             dum1 = dum1 * (mue(i).C-q_mue)^2;
%         end
%         dum2 = 0.0;
        minCFUE = inf;
        for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            fbs = fbs.SINR(SINR_FUE_Vec(j));
            if minCFUE>fbs.C_FUE
                minCFUE = fbs.C_FUE;
            end
            %dum2 = dum2 + (fbs.C_FUE-q_fue)^2;
            FBS{j}=fbs;
        end
        
%         R = K - dum2*dum1;
        for j=1:size(FBS,2)
            fbs = FBS{j};
            qMax=max(QTable,[],2);
            for jjj = 1:31
                if actions(1,jjj) == fbs.P
                    break;
                end
            end
            for kk = 1:size(states,1)
                if states(kk,:) == fbs.state
                    break;
                end
            end
            % CALCULATING NEXT STATE AND REWARD
            %beta = fbs.C_profile(episode)/fbs.C_FUE; %
            switch rewards
                case 'double_rate' 
                     if (fbs.C_FUE>q_fue)
                         R = 2*fbs.C_FUE;
                         fbs = fbs.reward(R);
        %                 R = beta*fbs.C_FUE;
                     else
                        R = fbs.C_FUE;
                        fbs = fbs.reward(R);
        %                 R = beta*fbs.C_FUE - (1/beta)*dum1*(minCFUE-q_fue).^2;
        %                 R = beta*fbs.C_FUE - (1/beta)*dum1;
                     end
                     
                case 'distance'
                    beta =fbs.d_near/dth; 
                    R= beta*fbs.C_FUE-(1/beta)*(fbs.C_FUE-q_fue)^2;
                    fbs = fbs.reward(R);
             end
%             if R<0
%                 R=0;
%             end
            for nextState=1:size(states,1)
                if states(nextState,:) == fbs.state
                    QTable(kk,jjj) = QTable(kk,jjj) + alpha*(R+gamma*qMax(nextState)-QTable(kk,jjj));
                end
            end
            FBS{j}=fbs;
        end
         %accumulated_reward(episode)=R;
%          if episode==1
%               accumulated_reward(episode)=R;
%          else
%              accumulated_reward(episode)=R+accumulated_reward(episode-1);
%          end
%          
        accumulated_reward(episode)=R;
        % break if convergence: small deviation on q for 1000 consecutive
        errorVector(episode) =  sum(sum(abs(Q1-QTable)));
        if sum(sum(abs(Q1-QTable)))<0.001 && sum(sum(QTable >0)) % Find out the deviation of Q-Table from the previous values, it means trying to find the optimal Q-Tables
            if count>1000
                episode  % report last episode
                break % for
            else
                count=count+1; % set counter if deviation of q is small
            end
        else
            Q1=QTable;
            count=0;  % reset counter when deviation of q from previous q is large
        end

%         if selectedMUE.C < gamma_th
%             I = 1;
%         else
%             I = 0;
%         end
% 
%         for j=1:size(FBS,2) 
%             fbs = FBS{j};
%             fbs.state(1,1) = I;
%             FBS{j} = fbs;
%         end
           case 'Exhaustive Search'
          textprogressbar((episode/Iterations)*100);
%          step_size = 46/size(FBS,2);
%         permutedPowers = -20:step_size:25;
%         
%         
%         for i=1:length(permutedPowers)
%             for j=1:length(permutedPowers)
%                 power(i,j)=permutedPowers(j);
%             end
%         end
%         horzon = horzcat(permutedPowers',power);
%         ee=[];
%         for jj=1:size(horzon,1)
%             cc= horzon(jj,:);
%             dd= unique(cc,'stable');
%             ee= [ee;dd];
%         end
           if count>=size(Combin,1)
              count=1;
          end
            % Action selection with epsilon=0.1
            %if rand<epsilon   % For implementation of epsilon greedy algorithm (It will explore with the probability of epsilon) : Exploration here means that randomly search for other power options
                for j=1:size(FBS,2)
                    fbs = FBS{j};
                    fbs = fbs.setPower(Combin(count,j));
                    FBS{j} = fbs;
                end
             count=count+1;
          

        % calc FUEs and MUEs capacity
        SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, noise,small_scale);
        C_FUE_Vec = log2(1+SINR_FUE_Vec);
%         for i=1:size(mue,2)
%             MUE = mue(i);
%             MUE.SINR = SINR_MUE_4(G, L, FBS, MBS, MUE, noise);
%             MUE = MUE.setCapacity(log2(1+MUE.SINR));
%             mue(i)=MUE;
%         end
        
%         MUE_C(1,episode) = selectedMUE.C;
        xx(1,episode) = episode;
%         R = K - (selectedMUE.SINR - sinr_th)^2;
%             deviation_FUE=0.0;
%             for i=1:size(FBS,2)
%                 deviation_FUE = deviation_FUE + (fbs.C_FUE-q_M)^2;
%             end
%         dum1 = 1.0;
%         for i=1:size(mue,2)
%             dum1 = dum1 * (mue(i).C-q_mue)^2;
%         end
%         dum2 = 0.0;
        minCFUE = inf;
        for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            fbs = fbs.SINR(SINR_FUE_Vec(j));
            if minCFUE>fbs.C_FUE
                minCFUE = fbs.C_FUE;
            end
            %dum2 = dum2 + (fbs.C_FUE-q_fue)^2;
            FBS{j}=fbs;
        end
        
             sumrate1 =0;
%             temp_sumrate= sumrate;
            for i=1:size(FBS,2)
                fbs = FBS{i};
                sumrate1 = sumrate1+ fbs.C_FUE;
            end
            temp_sumrate(episode)= sumrate1;
            if episode>1
                if sumrate1>temp_sumrate(episode-1)
                   sumrate= sumrate1;
                   powervect = Combin(episode,:); 
                else
                    sumrate=sumrate;
                     powervect = powervect; 
                end
            else
               sumrate =sumrate1;
                powervect = Combin(episode,:);
            end
%         R = K - dum2*dum1;
%         for j=1:size(FBS,2)
%             fbs = FBS{j};
%             qMax=max(QTable,[],2);
%             for jjj = 1:31
%                 if actions(1,jjj) == fbs.P
%                     break;
%                 end
%             end
%             for kk = 1:size(states,1)
%                 if states(kk,:) == fbs.state
%                     break;
%                 end
%             end
%             % CALCULATING NEXT STATE AND REWARD
%             %beta = fbs.C_profile(episode)/fbs.C_FUE; %
%             switch rewards
%                 case 'double_rate' 
%                      if (fbs.C_FUE>q_fue)
%                          R = 2*fbs.C_FUE;
%         %                 R = beta*fbs.C_FUE;
%                      else
%                         R = fbs.C_FUE;
%         %                 R = beta*fbs.C_FUE - (1/beta)*dum1*(minCFUE-q_fue).^2;
%         %                 R = beta*fbs.C_FUE - (1/beta)*dum1;
%                      end
%                      
%                 case 'distance'
%                     beta =fbs.d_near/dth; 
%                     R= beta*fbs.C_FUE-(1/beta)*(fbs.C_FUE-q_fue)^2;
%              end
% %             if R<0
% %                 R=0;
% %             end
%             for nextState=1:size(states,1)
%                 if states(nextState,:) == fbs.state
%                     QTable(kk,jjj) = QTable(kk,jjj) + alpha*(R+gamma*qMax(nextState)-QTable(kk,jjj));
%                 end
%             end
%             FBS{j}=fbs;
%         end
%          accumulated_reward(episode)=R;
%         % break if convergence: small deviation on q for 1000 consecutive
%         errorVector(episode) =  sum(sum(abs(Q1-QTable)));
%         if sum(sum(abs(Q1-QTable)))<0.001 && sum(sum(QTable >0)) % Find out the deviation of Q-Table from the previous values, it means trying to find the optimal Q-Tables
%             if count>1000
%                 episode  % report last episode
%                 break % for
%             else
%                 count=count+1; % set counter if deviation of q is small
%             end
%         else
%             Q1=QTable;
%             count=0;  % reset counter when deviation of q from previous q is large
%         end      
               
        case 'Waterfilling'
           
            %% waterfilling algorithm for power allocation
             %alocate waterfilling based power among the UAVs, that is,
            %allocate power based on their channel and  situation
            
            %Pn=(L(1:fbsCount,1)').^2;
%             P = 10.^((Pmax-30)/10);
%             P = waterfill(Pt,Pn);
%             
%             for j=1:size(FBS,2)
%                     fbs = FBS{j};
%                     fbs = fbs.setPower(10*log(P));
%                     FBS{j} = fbs;
%             end
            [SINR_FUE_Vec, rx_inter] = SINR_FUE_2(G, L, FBS, MBS, noise,small_scale);
            % Step to update the powers by including the effect of
            % interference as well. This will more realistically allocate
            % power among users
            %Pn=(G(1:fbsCount,1)'+rx_inter);
            for jj=1:size(FBS,2)
            %Pn=G(jj,jj)+rx_inter(1,jj);
            fbs = FBS{j};
            Pn=G(jj,jj)+rx_inter(1,jj);
            %Pt = 10.^((Pmax-30)/10);
            Pt = 10.^((Pmax-30)/10);
            P(jj) = waterfill(Pt,Pn);
            end
            for j=1:size(FBS,2)
                    fbs = FBS{j};
                    fbs = fbs.setPower(10*log10(P(j))+30);
                    FBS{j} = fbs;
            end
            
            minCFUE = inf;
        
            for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            if minCFUE>fbs.C_FUE
                minCFUE = fbs.C_FUE;
            end
            %dum2 = dum2 + (fbs.C_FUE-q_fue)^2;
            FBS{j}=fbs;
            end
        case 'Max_power'
            for j=1:size(FBS,2)
                    fbs = FBS{j};
                    fbs = fbs.setPower(Pmax);
                    FBS{j} = fbs;
            end
            [SINR_FUE_Vec, rx_inter] = SINR_FUE_2(G, L, FBS, MBS, noise,small_scale);
            minCFUE = inf;
        
            for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            if minCFUE>fbs.C_FUE
                minCFUE = fbs.C_FUE;
            end
            %dum2 = dum2 + (fbs.C_FUE-q_fue)^2;
            FBS{j}=fbs;
            end
            
        case 'PSO'
            textprogressbar((episode/Iterations)*100);
%               for j=1:size(FBS,2)
%                     fbs = FBS{j};
%                      positions=PSO(FBS,MBS,mue,NumRealization,fbsCount);
%                     fbs.X = x;
%                      fbs.Y = y;
%                     FBS{j} = fbs;
%             end
%              
%             
%               ii=1;
%               jj=2;
%               for a=1:size(FBS,2)
%                   fbs = FBS{a};
%                   fbs.X=positions(:,ii);
%                   fbs.Y=positions(:,jj);
%                     ii=ii+2;
%                     jj=jj+2;
%                   FBS{a}=fbs;
% 
%               end

           if count>=size( Powers,1)
              count=1;
          end
            % Action selection with epsilon=0.1
            %if rand<epsilon   % For implementation of epsilon greedy algorithm (It will explore with the probability of epsilon) : Exploration here means that randomly search for other power options
                for j=1:size(FBS,2)
                    fbs = FBS{j};
                    fbs = fbs.setPower(Powers(count,j));
                    FBS{j} = fbs;
                end
             count=count+1;
%            for j=1:size(FBS,2)
%                     fbs = FBS{j};
%                     fbs = fbs.setPower(Pmax);
%                     FBS{j} = fbs;
%             end

%         permutedPowers = randperm(size(actions,2),size(FBS,2));
%     
%     % y=randperm(size(permutedPowers,1));
%         for j=1:size(FBS,2)
%             fbs = FBS{j};
%             fbs = fbs.setPower(actions(permutedPowers(j)));
%             FBS{j} = fbs;
%         end
            %[G,L]= measure_channel_UAV_PSO1(FBS,MBS,mue,NumRealization, varying_threshold);
            [SINR_FUE_Vec, rx_inter]= SINR_FUE_2(G, L, FBS, MBS, noise,small_scale);
            %[SINR_FUE_Vec, rx_inter] = SINR_FUE_PSO(G, L, FBS, MBS, noise);
            minCFUE = inf;
        
            for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            if minCFUE>fbs.C_FUE
                minCFUE = fbs.C_FUE;
            end
            %dum2 = dum2 + (fbs.C_FUE-q_fue)^2;
            FBS{j}=fbs;
            end
            
            
            
           case 'PSO-PA'
             power=PSO_PA(FBS,MBS,mue,NumRealization,fbsCount,h,10,rad_Ubs,actions,Pmin,step,Pmax);
             textprogressbar((episode/Iterations)*100);
                
              for a=1:size(FBS,2)
                  fbs = FBS{a};
                  fbs= fbs.setPower(power(:,a));
                  fbs1 = FemtoStation_3S(fbs.X,fbs.Y, h, MBS, mue, 10,rad_Ubs);
                  fbs.FUEX = fbs1.FUEX;
                  fbs.FUEY = fbs1.FUEY;
                  FBS{a}=fbs;

              end
              [G,L]= measure_channel_UAV(FBS,MBS,mue,NumRealization, varying_threshold);
              
              [SINR_FUE_Vec, rx_inter]= SINR_FUE_2(G, L, FBS, MBS, noise,small_scale);
            %[SINR_FUE_Vec, rx_inter] = SINR_FUE_PSO(G, L, FBS, MBS, noise);
            minCFUE = inf;
        
            for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            if minCFUE>fbs.C_FUE
                minCFUE = fbs.C_FUE;
            end
            %dum2 = dum2 + (fbs.C_FUE-q_fue)^2;
            FBS{j}=fbs;
            end
            
       case'EPA'
           
            textprogressbar((episode/Iterations)*100);
            
%             for i=1:length(actions)
%                  Powers(i,:) = repelem(actions(i),fbsCount);
%             end
           if count>=size( Powers,1)
              count=1;
          end
            % Action selection with epsilon=0.1
            %if rand<epsilon   % For implementation of epsilon greedy algorithm (It will explore with the probability of epsilon) : Exploration here means that randomly search for other power options
                for j=1:size(FBS,2)
                    fbs = FBS{j};
                    fbs = fbs.setPower(Powers(count,j));
                    FBS{j} = fbs;
                end
             count=count+1;
          

        % calc FUEs and MUEs capacity
        SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, noise,small_scale);
        C_FUE_Vec = log2(1+SINR_FUE_Vec);
%         for i=1:size(mue,2)
%             MUE = mue(i);
%             MUE.SINR = SINR_MUE_4(G, L, FBS, MBS, MUE, noise);
%             MUE = MUE.setCapacity(log2(1+MUE.SINR));
%             mue(i)=MUE;
%         end
        
%         MUE_C(1,episode) = selectedMUE.C;
        xx(1,episode) = episode;
%         R = K - (selectedMUE.SINR - sinr_th)^2;
%             deviation_FUE=0.0;
%             for i=1:size(FBS,2)
%                 deviation_FUE = deviation_FUE + (fbs.C_FUE-q_M)^2;
%             end
%         dum1 = 1.0;
%         for i=1:size(mue,2)
%             dum1 = dum1 * (mue(i).C-q_mue)^2;
%         end
%         dum2 = 0.0;
        minCFUE = inf;
        for j=1:size(FBS,2)
            fbs = FBS{j};
            fbs = fbs.setCapacity(log2(1+SINR_FUE_Vec(j)));
            fbs = fbs.SINR(SINR_FUE_Vec(j));
            if minCFUE>fbs.C_FUE
                minCFUE = fbs.C_FUE;
            end
            %dum2 = dum2 + (fbs.C_FUE-q_fue)^2;
            FBS{j}=fbs;
        end
        
             sumrate1 =0;
%             temp_sumrate= sumrate;
            for i=1:size(FBS,2)
                fbs = FBS{i};
                sumrate1 = sumrate1+ fbs.C_FUE;
            end
            temp_sumrate(episode)= sumrate1;
            
            if episode>1
                if sumrate1>temp_sumrate(episode-1)
                   sumrate= sumrate1;
                   powervect = Powers(count,:);
                   powervector(episode,:)= Powers(count,:);
                else
                    sumrate=sumrate;
                     powervect = powervect;
                     powervector(episode,:)= Powers(count,:);
                end
            else
               sumrate =sumrate1;
                powervect = Powers(count,:);
                powervector(episode,:)= Powers(count,:);
            end
           
           
           
           
          
        end
    end
    [max_rate,pos]=max(temp_sumrate,[],2);
    power=powervector(pos,:);
      sum_R = zeros(1,fbsCount);  
    switch power_allocation_algo
        case 'EPA'
            min_CFUE = inf;
                for j=1:size(FBS,2)
                    C = FBS{1,j}.C_profile;
                   H=FBS{1,j}.h;
                    c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));   
                    h_(1,j)=H;
                    if min_CFUE > c_fue(1,j)
                        min_CFUE = c_fue(1,j);
                    end
                end
                answer.sum_CFUE = max_rate; 
                answer.powervectf=power;
                answer.min_CFUE = min_CFUE;
                answer.FBS = FBS;
                answer.episode = episode;
                answer.time = toc;
                Q = QTable;
                answer.Q = QTable;
                QFinal = answer;
            
        case 'Exhaustive Search'
            
            
                min_CFUE = inf;
                for j=1:size(FBS,2)
                    C = FBS{1,j}.C_profile;
                    H=FBS{1,j}.h;
                    c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));   
                    h_(1,j)=H;
                    if min_CFUE > c_fue(1,j)
                        min_CFUE = c_fue(1,j);
                    end
                end
                answer.sum_CFUE = sumrate; 
                answer.powervectf=powervect;
                answer.min_CFUE = min_CFUE;
                answer.FBS = FBS;
                answer.episode = episode;
                answer.time = toc;
                Q = QTable;
                answer.Q = QTable;
                QFinal = answer;
               
                 
        otherwise
      
                Q = QTable;
                answer.mue = mue;
                answer.Q = QTable;
                answer.Error = errorVector;
                answer.FBS = FBS;
                min_CFUE = inf;
                for j=1:size(FBS,2)
                    C = FBS{1,j}.C_profile;
                     sum_R(1,j)=sum(FBS{1,j}.R_profile);
                    H=FBS{1,j}.h;
                    switch power_allocation_algo
                        case 'Q_learning'   
                        %c_fue(1,j) = sum(C(40000:size(C,2)))/(-40000+size(C,2));
                        c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));
                        case 'Waterfilling'
                        c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));
                        case 'Exhaustive Search'
                        c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));
                        case 'Max_power'
                        c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));
                        case 'PSO'
                        c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));
                        case 'PSO-PA'
                        c_fue(1,j) = sum(C(1:size(C,2)))/(-1+size(C,2));
                    end
                    
                    h_(1,j)=H;
                    if min_CFUE > c_fue(1,j)
                        min_CFUE = c_fue(1,j);
                    end
                end
                sum_CFUE = 0.0;
                for i=1:size(FBS,2)
                    sum_CFUE = sum_CFUE + c_fue(1,i);
                end
                answer.AvgR = sum(sum_R)/fbsCount;
                answer.C_FUE = c_fue;
                answer.h_com = h_;
                answer.sum_CFUE = sum_CFUE;
                answer.min_CFUE = min_CFUE;
                answer.episode = episode;
                answer.time = toc;
                QFinal = answer;
    end
    %save(sprintf('results/pro_%d_%d_%d.mat',fbsCount, saveNum, h1),'QFinal');
end




