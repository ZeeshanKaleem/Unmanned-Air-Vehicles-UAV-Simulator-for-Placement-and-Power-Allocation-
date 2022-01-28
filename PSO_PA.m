function power = PSO_PA(fbs,MBS,mue,NumRealization,fbscount,h,cc,rad_Ubs,actions,Pmin,step,Pmax)
% tic 
% clc 
% clear all 
% close all
% rng default
varying_threshold=0;
fbsc=fbscount;
noise=-174;
small_scale=0;
range = 500;

% LB=[0 0;0 0;0 0; 0 0];% here each row represents each ABS (x,y) lower bound positions. 
% UB = [500 500;500 500;500 500; 500 500]; % Case of 4 ABS
if fbscount==4
LB=[-20;-20;-20;-20];% here each row represents each ABS (x,y) lower bound positions. 
UB = [25;25;25;25]; % Case of 4 ABS

elseif fbscount==8
LB=[-20;-20;-20;-20;-20;-20;-20;-20];% here each row represents each ABS (x,y) lower bound positions. 
UB = [25;25;25;25;25;25;25;25]; % Case of 8 ABS

elseif fbscount==16
    

end
       
%LB=[0 0 0];
%UB=[10 10 10];      %upper bounds of variables
% pso parameters values
m=1;            % number of variables
n=50;          % population size
wmax=0.8;       % inertia weight
wmin=0.6;       % inertia weight 
c1=2;           % acceleration factor
c2=2;           % acceleration factor
%X0=cell(1,fbsc);
X0=[];
fbsXx =[];
fbsYy =[];
SINR =[];
jointfbsXY =[];
Snr= zeros(n,fbscount);
cap=zeros(n,fbscount);
x=cell(1,fbsc);
v=cell(1,fbsc);
pbest =cell(1,fbsc);
gbest =cell(1,fbsc);


% pso main program----------------------------------------------------start 
maxite=100;    % set maximum number of iteration 
maxrun=1;      % set maximum number of runs need to be
ffmin = cell(1,fbsc);
ffite = cell(1,fbsc);
rgbest=cell(1,maxrun);
fff=cell(1,maxrun);
fbsX=cell(1,fbsc);
fbsY=cell(1,fbsc);
for run=1:maxrun
    %run
  for a=1:fbscount
    for i=1:n
        for j=1:m
            x0(i,j)=round(LB(a,j)+rand()*(UB(a,j)-LB(a,j)));
        end
    end
    X0=[X0 x0];
  end
    

    for a=1:fbscount
       
        fbsX{a}=X0(:,a);
           
    end
  for a=1:fbscount
      fbsXx= [fbsXx fbsX{a}];
    
  end
   for a=1:fbscount
      fbsXx1= fbsX{a};
      %fbsYy1= fbsY{a};
      jointfbsXY1 = fbsXx1;
      x{a}=jointfbsXY1; % initial population
      v{a}=0.1*jointfbsXY1; % initial velocity
      %v= [v vv]; % initial velocity
      jointfbsXY = jointfbsXY1;
  end
  
      
      for i=1:n
           for a=1:fbscount
               fbs{a}.P=fbsXx(i,a);
               %fbs{a}.Y=fbsYy(i,a);
               
                fbs1 = FemtoStation_3S(fbs{a}.X,fbs{a}.Y, h, MBS, mue, cc,rad_Ubs);
                fbs{a}.FUEX = fbs1.FUEX;
                fbs{a}.FUEY = fbs1.FUEY;
                
           end
           [G,L]= measure_channel_UAV_PSO1(fbs,MBS,mue,NumRealization, varying_threshold);
           Snr(i,:)= SINR_FUE_PSO(G, L, fbs, MBS, noise,small_scale);
           f0(i,:) = ofun1(Snr(i,:));
      end
  
  
%   x=x0;       % initial population
%     v=0.1*x0;   % initial velocity 
%     for i=1:n
%         f0(i,1)=ofun(x0(i,:),fbs,MBS,mue,NumRealization);
%     end
    [fmin0,index0]=min(f0);
    for a=1:fbscount
    pbest1=fbsXx(:,a);               % initial pbest
    %pbest2=fbsYy(:,a);
    pbest{a}= pbest1;
    gbest{a}=pbest{a}(index0(a),:);     % initial gbest
    end
    % pso initialization------------------------------------------------end
% pso algorithm---------------------------------------------------start
ite=1;
tolerance=1;
while ite<=maxite && tolerance>10^-12
w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight
% pso velocity updates
for a=1:fbscount
for i=1:n
for j=1:m
v{a}(i,j)=w*v{a}(i,j)+c1*rand()*(pbest{a}(i,j)-x{a}(i,j))...
+c2*rand()*(gbest{a}(1,j)-x{a}(i,j));
end
end
end
% pso position update
for a=1:fbscount
for i=1:n
for j=1:m
x{a}(i,j)=x{a}(i,j)+v{a}(i,j);
end
end
end
% handling boundary violations
for a=1:fbscount
for i=1:n
for j=1:m
if x{a}(i,j)<LB(j)
x{a}(i,j)=LB(j);
elseif x{a}(i,j)>UB(j)
x{a}(i,j)=UB(j);
end
end
end
end
% evaluating fitness

  for a=1:fbscount
       
        fbsX{a}=x{a}(:,1);
        %fbsY{a}=x{a}(:,2);
           
  end


     for i=1:n
           for a=1:fbscount
               fbs{a}.P=fbsX{a}(i,1);
               fbs1 = FemtoStation_3S(fbs{a}.X,fbs{a}.Y, h, MBS, mue, cc,rad_Ubs);
                fbs{a}.FUEX = fbs1.FUEX;
                fbs{a}.FUEY = fbs1.FUEY;
           end
           [G,L]= measure_channel_UAV_PSO1(fbs,MBS,mue,NumRealization, varying_threshold);
           Snr1(i,:)= SINR_FUE_PSO(G, L, fbs, MBS, noise,small_scale);
           f(i,:) = ofun1(Snr1(i,:));
      end


% for i=1:n
% f(i,1)=ofun(x(i,:));
% end
% updating pbest and fitness
for a=1:fbscount
for i=1:n
if f(i,1)<f0(i,1)
    pbest{a}(i,:)=x{a}(i,:);
f0(i,1)=f(i,1);
end
end
end
[fmin,index]=min(f0); % finding out the best particle
for a=1:fbscount
ffmin{a}(ite,run)=fmin(a); % storing best fitness
ffite{a}(run)=ite; % storing iteration count
end
% updating gbest and best fitness
for a=1:fbscount
if fmin(a)<fmin0(a)
gbest{a}=pbest{a}(index(a),:);
fmin0(a)=fmin(a);
end
end
% calculating tolerance
for a=1:fbscount
if ite>100
tolerance=abs(ffmin{a}(ite-100,run)-fmin0(a));
end
% displaying iterative results
if ite==1
%disp(sprintf('Iteration Best particle Objective fun'));
end
%disp(sprintf('%8g %8g %8.4f',ite,index,fmin0(a)));
end
ite=ite+1;
end
% pso algorithm-----------------------------------------------------end
gbest;
%fvalue=10*(gbest(1)-1)^2+20*(gbest(2)-2)^2+30*(gbest(3)-3)^2;

   for a=1:fbscount
       FBS = fbs{a};
       FBS.P=gbest{a}(:,1);
       %FBS.Y=gbest{a}(:,2);
       fbs1 = FemtoStation_3S(FBS.X,FBS.Y, h, MBS, mue, cc,rad_Ubs);
                FBS.FUEX = fbs1.FUEX;
                FBS.FUEY = fbs1.FUEY;
       fbs{a}=FBS;
       
   end
   [G,L]= measure_channel_UAV_PSO1(fbs,MBS,mue,NumRealization, varying_threshold);
   Snr2 = SINR_FUE_PSO(G, L, fbs, MBS, noise,small_scale);
   fvalue = ofun1(Snr2);
   sum_rate = sum(fvalue);

fff{:,run}=sum_rate;
rgbest{:,run}=gbest;
%disp(sprintf('--------------------------------------'));
end
% pso main program------------------------------------------------------end
%disp(sprintf('\n'));
%disp(sprintf('*********************************************************'));
%disp(sprintf('Final Results-----------------------------'));
[bestfun,bestrun]=min([fff{:,:}]);
best_variables=rgbest{:,bestrun};
%disp(sprintf('*********************************************************'));
power= [best_variables{:,:}];
%toc

% PSO convergence characteristic
% plot(ffmin{1}(1:ffite{1,1}(bestrun),bestrun),'-k');
% xlabel('Iteration');
% ylabel('Fitness function value');
% title('PSO convergence characteristic')
%##########################################################################