function [xUAV,yUAV] = optimal_x_y_points(num_points,range,fbscount, ploting)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if fbscount==4
x= [randi([1 range],num_points,1); randi([-range 1],num_points,1); randi([-range 1],num_points,1); randi([1 range],num_points,1)];
y =[randi([1 range],num_points,1);  randi([1 range],num_points,1); randi([-range 1],num_points,1); randi([-range 1],num_points,1)];
elseif fbscount==8
    range = range/2;  
    x= [randi([1 range],num_points,1); randi([range range+range],num_points,1); randi([-range 1],num_points,1); randi([-range-range -range],num_points,1); ...
        randi([-range 1],num_points,1); randi([-range-range -range],num_points,1); randi([1 range],num_points,1);randi([range range+range],num_points,1)];
    
    y =[randi([1 range+range],num_points,1);randi([1 range+range],num_points,1);  randi([1 range+range],num_points,1); randi([1 range+range],num_points,1); ...
        randi([-range-range 1],num_points,1);randi([-range-range 1],num_points,1); randi([-range-range 1],num_points,1); randi([-range-range 1],num_points,1)];
elseif fbscount==16
    range = range/4;
    r=range;% 25
    r1=r+range; % 50
    r2 = r1+range; % 75
    r3 = r2+range; % 100
    
    
    x= [randi([1 r],num_points,1); randi([r r1],num_points,1);randi([r1 r2],num_points,1); randi([r2 r3],num_points,1);...
        randi([-r 1],num_points,1); randi([-r1 -r],num_points,1);randi([-r2 -r1],num_points,1); randi([-r3 -r2],num_points,1);...
        randi([-r 1],num_points,1); randi([-r1 -r],num_points,1);randi([-r2 -r1],num_points,1); randi([-r3 -r2],num_points,1);...
        randi([1 r],num_points,1); randi([r r1],num_points,1);randi([r1 r2],num_points,1); randi([r2 r3],num_points,1)];
     
    y =[randi([1 4*r],num_points,1);randi([1 4*r],num_points,1);  randi([1 4*r],num_points,1); randi([1 4*r],num_points,1); ...
        randi([1 4*r],num_points,1);randi([1 4*r],num_points,1);  randi([1 4*r],num_points,1); randi([1 4*r],num_points,1);...
        randi([-4*r 1],num_points,1);randi([-4*r 1],num_points,1);  randi([-4*r 1],num_points,1); randi([-4*r 1],num_points,1);...
        randi([-4*r 1],num_points,1);randi([-4*r 1],num_points,1);  randi([-4*r 1],num_points,1); randi([-4*r 1],num_points,1)];
else 
    disp('not supported')
end

X = [x y];
opts = statset('Display','final');
[cidx, ctrs] = kmeans(X, fbscount, 'Distance','city', ...[cidx, ctrs] = kmeans(X, 2, 'Distance','city', ...
                              'Replicates',5, 'Options',opts);

d=ctrs';
xUAV = d(1,:);
yUAV = d(2,:);
                          
if ploting                         
figure
if fbscount==4
plot(X(cidx==1,1),X(cidx==1,2),'r.')
hold on
plot(X(cidx==2,1),X(cidx==2,2),'b.')
hold on
plot(X(cidx==3,1),X(cidx==3,2),'g*')
hold on
plot(X(cidx==4,1),X(cidx==4,2),'y*')
plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize' ,15, 'LineWidth' ,3)
elseif fbscount==8
    plot(X(cidx==1,1),X(cidx==1,2),'r.')
hold on
plot(X(cidx==2,1),X(cidx==2,2),'b.')
hold on
plot(X(cidx==3,1),X(cidx==3,2),'g*')
hold on
plot(X(cidx==4,1),X(cidx==4,2),'y*')
plot(X(cidx==5,1),X(cidx==5,2),'k*')
plot(X(cidx==6,1),X(cidx==6,2),'m*')
plot(X(cidx==7,1),X(cidx==7,2),'c*')
plot(X(cidx==8,1),X(cidx==8,2), 'r*')
plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize' ,15, 'LineWidth' ,3)
elseif fbscount==16
plot(X(cidx==1,1),X(cidx==1,2),'r.')
hold on
plot(X(cidx==2,1),X(cidx==2,2),'b.')
hold on
plot(X(cidx==3,1),X(cidx==3,2),'g*')
hold on
plot(X(cidx==4,1),X(cidx==4,2),'y*')
plot(X(cidx==5,1),X(cidx==5,2),'k*')
plot(X(cidx==6,1),X(cidx==6,2),'m*')
plot(X(cidx==7,1),X(cidx==7,2),'c*')
plot(X(cidx==8,1),X(cidx==8,2), 'r*')
plot(X(cidx==9,1),X(cidx==9,2), 'r*')
plot(X(cidx==10,1),X(cidx==10,2), 'color', [.5 .6 .7], 'marker','.')
plot(X(cidx==10,1),X(cidx==10,2), 'color', [.4 .5 .7], 'marker','.')
plot(X(cidx==11,1),X(cidx==11,2), 'color', [.6 .7 .7], 'marker','.')
plot(X(cidx==12,1),X(cidx==12,2), 'color', [.5 .8 .7], 'marker','.')
plot(X(cidx==13,1),X(cidx==13,2), 'color', [.5 .6 .7], 'marker','*')
plot(X(cidx==14,1),X(cidx==14,2), 'color', [.3 .6 .7], 'marker','*')
plot(X(cidx==15,1),X(cidx==15,2), 'color', [.4 .6 .7], 'marker','*')
plot(X(cidx==16,1),X(cidx==16,2), 'color', [.7 .6 .7], 'marker','*')
plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize' ,15, 'LineWidth' ,3)
end
end
end

