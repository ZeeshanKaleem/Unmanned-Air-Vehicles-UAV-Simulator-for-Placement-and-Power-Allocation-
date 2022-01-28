function [G, L] = measure_channel_UAV(FBS,MBS,MUE,NumRealization, varying_threshold)
    fbsNum = size(FBS,2);
   if size(MUE,2)<2
        G = zeros(fbsNum+1, fbsNum+1);
        L = zeros(fbsNum+1, fbsNum+1);
   else 
        G = zeros(fbsNum+1, fbsNum+size(MUE,2)); % Matrix Containing small scale fading coefficients
        L = zeros(fbsNum+1, fbsNum+size(MUE,2)); % Matrix Containing large scale fading coefficients
   end
    %P_LOS = zeros(fbsNum+1, fbsNum+1);
    %pl_UAVUE_UAV = zeros(NUE,fbsNum);
%     NUE=1;
% %% old curves
%     f_c=0.9*10^9;
%     s_light= 3*10^8;
%     
% %     f=0.9; % Frequency in GHz (It means they took 900MHz freq.)
% %     PLi = -1.8*f^2+10.6*f-5.5; % This value is taken for Zerro N_w according to reference [22]
%     for i=1:fbsNum
%         xAgent = FBS{i}.X;
%         yAgent = FBS{i}.Y;
%         hAgent = FBS{i}.h;
%         A=FBS{i}.eta_LOS-FBS{i}.eta_NLOS;
%         B=20*(log((4*pi*f_c)/s_light))+FBS{i}.eta_NLOS;
%         a=FBS{i}.a;
%         b=FBS{i}.b;
%         % %% LOS probability calculations
%         %P_LOS = zeros(NUE,fbsNum); % LOS probability
% 
%         for j=1:fbsNum
%             d = sqrt((xAgent-FBS{j}.FUEX)^2+(yAgent-FBS{j}.FUEY)^2);
%            
%                 %PL0 = 62.3+40*log10(d/5);
%                 if varying_threshold
%                     P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
%                     PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
%                 else
%                     P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
%                     PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
%                 end
%              L(i,j) = 10^((PL0)/10);  
%         end   
%            %%
%            %Plot height vs. pathloss
%            height=0:1:1200;
%            PL_h=zeros(1,length(height));
%            for h=1:length(height)
%                P_LOS=1./(1+(a.*exp(-b.*(atan( height(h)./d)-a))));
%                PL_h(:,h)=20*(log(sqrt(height(h)^2+d^2)/cos(atan( height(h)./d))))+A*P_LOS+B;
%            
%            end
%            plot(height,PL_h)
           
    %% new curves
    
    f_c=0.9e+9;
    s_light= 3*10^8;
    
%     f=0.9; % Frequency in GHz (It means they took 900MHz freq.)
%     PLi = -1.8*f^2+10.6*f-5.5; % This value is taken for Zerro N_w according to reference [22]
    for i=1:fbsNum
        xAgent = FBS{i}.X;
        yAgent = FBS{i}.Y;
        hAgent = FBS{i}.h;
        A=FBS{i}.eta_LOS-FBS{i}.eta_NLOS;
        B=20*(log((4*pi*f_c)/s_light))+FBS{i}.eta_NLOS;
        a=FBS{i}.a;
        b=FBS{i}.b;
        % %% LOS probability calculations
        %P_LOS = zeros(NUE,fbsNum); % LOS probability

        for j=1:fbsNum
            d = sqrt((xAgent-FBS{j}.FUEX)^2+(yAgent-FBS{j}.FUEY)^2);
           
                %PL0 = 62.3+40*log10(d/5);
                if varying_threshold
                    P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
                    PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
                else
                    theta =atan(hAgent/d);
               denom = 1+(a*exp(-b*(((180/pi)*theta)-a)));
               P_LOS = 1/denom;
               hei= sqrt((hAgent)^2+(d)^2);
%                tem=4*pi*f_c*hei/s_light;
%                PL_h(:,h)=20*log(tem)+P_LOS*eta_LOS+((1-P_LOS)*eta_NLOS);
                tem=4*pi*f_c/s_light;
               PL0=20*log(tem)+20*log(hei)+P_LOS*FBS{i}.eta_LOS+((1-P_LOS)*FBS{i}.eta_NLOS);
                    
                    
%                     P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
%                     PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
                end
             L(i,j) = 10^((PL0)/10);  
        end   
           %%
           %Plot height vs. pathloss
%            height=0:1:1200;
%            PL_h=zeros(1,length(height));
%            for h=1:length(height)
%                P_LOS=1./(1+(a.*exp(-b.*(atan( height(h)./d)-a))));
%                PL_h(:,h)=20*(log(sqrt(height(h)^2+d^2)/cos(atan( height(h)./d))))+A*P_LOS+B;
%            
%            end
%            plot(height,PL_h)
    
                   
           
           
           if size(MUE,2)>1
               for k=1:size(MUE,2)
               d = sqrt((xAgent-MUE(k).X)^2+(yAgent-MUE(k).Y)^2);
               %d = nearest_MUE(xAgent, yAgent, MUE);  
                    if varying_threshold
                        P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
                        PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
                    else
                        P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
                        PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
                    end
                        L(i,fbsNum+k) = 10^((PL0)/10); 
               end
                 
           else
               d = sqrt((xAgent-MUE.X)^2+(yAgent-MUE.Y)^2);
           
                if varying_threshold
                    P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
                    PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
                else
                    P_LOS=1./(1+(a.*exp(-b.*(180/pi*atan( hAgent./d)-a))));
                    PL0=20*(log(sqrt(hAgent^2+d^2)))+A*P_LOS+B;
                end
                        
                 L(i,fbsNum+1) = 10^((PL0)/10);        
            end
        
        
%         PL0 = 62.3+32.*log10(d/5);
%         L(i,fbsNum+1) = 10.^((PL0 + PLi)/10);
%         
        d = sqrt((MBS.X-FBS{i}.FUEX)^2+(MBS.Y-FBS{i}.FUEY)^2);
        PL_BS = 15.3+37.6*log10(d);
        L(fbsNum+1,i) = 10^((PL_BS)/10);
    end
    
 if size(MUE,2)<2
    d = sqrt((MBS.X-MUE.X).^2+(MBS.Y-MUE.Y).^2);
   % PL_BS = 15.3+37.6*log10(d);
    PL_BS= 34.5 + 35*log10(d);
    L(fbsNum+1,fbsNum+1) = 10.^((PL_BS)/10);
 else
     for kk=1:size(MUE,2)
         d = sqrt((MBS.X-MUE(kk).X).^2+(MBS.Y-MUE(kk).Y).^2);
   % PL_BS = 15.3+37.6*log10(d);
         PL_BS= 34.5 + 35*log10(d);
         L(fbsNum+1,fbsNum+kk) = 10.^((PL_BS)/10);
     end
 end
 
 if size(MUE,2)<2
    Hij = abs((1/sqrt(2)) * (randn(fbsNum+1, fbsNum+1, NumRealization)+1i*randn(fbsNum+1, fbsNum+1, NumRealization)));
    hij = Hij.^2;
    for i=1:fbsNum+1
        for j=1:fbsNum+1
            G(i,j)=(sum(hij(i,j,:))/NumRealization);
        end
    end
 else
      Hij = abs((1/sqrt(2)) * (randn(fbsNum+1, fbsNum+size(MUE,2), NumRealization)+1i*randn(fbsNum+1, fbsNum+size(MUE,2), NumRealization)));
    hij = Hij.^2;
    for i=1:fbsNum+1
        for j=1:fbsNum+size(MUE,2)
            G(i,j)=(sum(hij(i,j,:))/NumRealization);
        end
    end
 end
     
end





% %% Distance calculation 
% %MUE to MBS horizontal distance 
% D_MUE_MBS = zeros(nMUE,nMBS); % home mue to home mbs
% x=0;
% for i = 1:nMBS
%     c=(xMUE(1+x:nMUE+x)-xMBS(i)).^2+(yMUE(1+x:nMUE+x)-yMBS(i)).^2;
%     x=x+nMUE;
%     D_MUE_MBS(:,i) = sqrt(c); 
% end
% 
% D_all_MUE_MBS = zeros(nMUE*nMBS,nMBS); % all mue to all mbs
% x=0;
% for i = 1:nMBS
%     c=(xMUE-xMBS(i)).^2+(yMUE-yMBS(i)).^2;
%     x=x+nMUE;
%     D_all_MUE_MBS(:,i) = sqrt(c); 
% end
% 
% %UAV-UE to UAV horizontal distance 
% D_UAVUE_UAV = zeros(NUvUE,nUAV); % home mue to home mbs
% x=0;
% for i = 1:nUAV
%     c=(xUvUE(1+x:NUvUE+x)-xUAV(i)).^2+(yUvUE(1+x:NUvUE+x)-yUAV(i)).^2;
% %     c=(xUvUE(1+x:NUvUE+x)-xUAV(:)).^2+(yUvUE(1+x:NUvUE+x)-yUAV(:)).^2;
%     x=x+NUvUE;
%     D_UAVUE_UAV(:,i) = sqrt(c); 
% end
% 
% %UAV-UE to MBS horizontal distance 
% % D_UAVUE_MBS = zeros(NUvUE,nMBS); % UAV UE to MBs
% % x=0;
% % for i = 1:nMBS
% %     c=(xUvUE(1+x:NUvUE+x)-xMBS(i)).^2+(yUvUE(1+x:NUvUE+x)-yMBS(i)).^2;
% %     x=x+NUvUE;
% %     D_UAVUE_MBS(:,i) = sqrt(c); 
% % end
% 
% % UAV-UE to all MBs
% D_all_UAVUE_mbs = zeros(NUvUE*nUAV*nMBS,nMBS); % all mue to all mbs
% for i = 1:nMBS
%     c=(xUvUE-xMBS(i)).^2+(yUvUE-yMBS(i)).^2;
%     D_all_UAVUE_mbs(:,i) = sqrt(c); 
% end
% 
% %MUE to UAV horizontal distance 
% D_MUE_UAV = zeros(nMUE,nUAV); % home mue to UAVs
% x=0;
% for i = 1:nUAV
%     c=(xMUE(1+x:nMUE+x)-xUAV(i)).^2+(yMUE(1+x:nMUE+x)-yUAV(i)).^2;
%     x=0;
%     D_MUE_UAV(:,i) = sqrt(c); 
% end
% 
% % MUE-all UAV
% D_all_mue_UAV = zeros(nMUE*nMBS,nUAV*nMBS); % all mue to all UAVs
% for i = 1:nUAV*nMBS
%     c=(xMUE-xUAV(i)).^2+(yMUE-yUAV(i)).^2;
%     D_all_mue_UAV(:,i) = sqrt(c); 
% end
% 
% % UAV-UE -all UAV
%     
% D_all_UAVUE_UAV = zeros(NUvUE*nUAV*nMBS,nUAV*nMBS); % all mue to all mbs
% for i = 1:nUAV*nMBS
%     c=(xUvUE-xUAV(i)).^2+(yUvUE-yUAV(i)).^2; % Each coloumn shows that all users in a system connects with the UAV1..so on
%     D_all_UAVUE_UAV(:,i) = sqrt(c); 
% end
% 
% 
% D1=[D_all_UAVUE_UAV D_all_UAVUE_mbs];
% D2=[D_all_mue_UAV D_all_MUE_MBS];
% D=vertcat(D1,D2);    % ditance between all users and all BSs  fifth coloumn is for macro BS
% 
% 
% %% LOS probability calculations
% 
% %UAV-UE to UAV horizontal distance 
% P_LOS = zeros(NUvUE,nUAV); % LOS probability
% for i = 1:nUAV
%     if varying_threshold
%         P_LOS(:,i)=1./(1+(a.*exp(-b.*(atan(h./D_UAVUE_UAV(:,i))-a))));
%     else
%         P_LOS(:,i)=1./(1+(a.*exp(-b.*(atan(h(l,:)./D_UAVUE_UAV(:,i))-a))));
%     end
% end
% 
% %% UAV-UE to UAV horizontal distance with all UAVs
% P_LOS_all = zeros(NUvUE*nUAV*nMBS,nUAV*nMBS); % LOS probability
% for i = 1:nUAV
%     if varying_threshold
%         P_LOS_all(:,i)=1./(1+(a.*exp(-b.*(atan(h./D_all_UAVUE_UAV(:,i))-a))));
%     else
%         P_LOS_all(:,i)=1./(1+(a.*exp(-b.*(atan(h(l,:)./D_all_UAVUE_UAV(:,i))-a))));
%     end
% end
% 
% %% MUE to UAV 
% P_LOS_all_mue = zeros(nMUE*nMBS,nUAV*nMBS); % LOS probability
% for i = 1:nUAV
%     if varying_threshold
%         P_LOS_all_mue(:,i)=1./(1+(a.*exp(-b.*(atan(h./D_all_mue_UAV(:,i)-a)))));
%     else
%         P_LOS_all_mue(:,i)=1./(1+(a.*exp(-b.*(atan(h(l,:)./D_all_mue_UAV(:,i)-a)))));
%     end
% end
% 
% %% Pathloss calculation for UAV-UE with UAVs
% pl_UAVUE_UAV = zeros(NUvUE,nUAV);
% A=eta_LOS-eta_NLOS;
% B=20*(log((4*pi*f_c)/s_light))+eta_NLOS;
% for i = 1:nUAV
%     if varying_threshold
%        pl_UAVUE_UAV(:,i)=20*(log(sqrt(h^2+(D_UAVUE_UAV(:,i)).^2)))+A*P_LOS(:,i)+B;
%     else
%         pl_UAVUE_UAV(:,i)=20*(log(sqrt(h(l,:)^2+(D_UAVUE_UAV(:,i)).^2)))+A*P_LOS(:,i)+B;
%     end
% end
% %% Pathloss calculation for UAV-UE with all UAVs
% pl_UAVUE_UAV_all = zeros(NUvUE*nUAV*nMBS,nUAV*nMBS);
% % A=eta_LOS-eta_NLOS;
% % B=20*log((4*pi*f_c)/s_light)+eta_NLOS;
% for i = 1:nUAV
%     if varying_threshold
%         pl_UAVUE_UAV_all(:,i)=20*(log(sqrt(h^2+D_all_UAVUE_UAV(:,i).^2)))+A*P_LOS_all(:,i)+B;
%     else
%         pl_UAVUE_UAV_all(:,i)=20*(log(sqrt(h(l,:)^2+D_all_UAVUE_UAV(:,i).^2)))+A*P_LOS_all(:,i)+B;
%     end
% end
% 
% %% Pathloss calculation for MUE with all UAVs
% pl_MUE_UAV_all = zeros(nMUE*nMBS,nUAV*nMBS);
% % A=eta_LOS-eta_NLOS;
% % B=20*log((4*pi*f_c)/s_light)+eta_NLOS;
% for i = 1:nUAV
%     if  varying_threshold
%         pl_MUE_UAV_all(:,i)=20*(log(sqrt(h^2+D_all_mue_UAV(:,i).^2)))+A*P_LOS_all_mue(:,i)+B;
%     else
%         pl_MUE_UAV_all(:,i)=20*(log(sqrt(h(l,:)^2+D_all_mue_UAV(:,i).^2)))+A*P_LOS_all_mue(:,i)+B;
%     end
% end
% 
% %% Pl calculation for MUE with MBs
% 
% pl_mue_mbs=15.3+37.6*log10(D_MUE_MBS);
% 
% %% Pl calculation for UAV-UE with MBs
% 
% pl_UAVUE_mbs=15.3+37.6*log10(D_all_UAVUE_mbs);
% 
% %% Total Pathloss calculation
% 
% pl_UAVUE_UAV_MBS= [pl_UAVUE_UAV_all pl_UAVUE_mbs];
% pl_MUE_UAV_MBS= [pl_MUE_UAV_all pl_mue_mbs];
% 
% path_loss = [pl_UAVUE_UAV_MBS; pl_MUE_UAV_MBS];