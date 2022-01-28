function SINR = SINR_FUE_PSO(G, L, FBS, MBS, sigma2,small_scale)
    MBS_P = 10^((MBS.P-30)/10); %Standarf formula for dBm to Watt conversion
    fbsNum = size(FBS,2);
    SINR = zeros(1,fbsNum);
    rx_inter = zeros(1,fbsNum);
    sigma = 10^((sigma2-30)/10);
    P_interf = 0.0;
    pAgent = zeros(1,fbsNum);
    for i=1:fbsNum
        pAgent(i) = 10.^((FBS{i}.P-30)/10);
        %pAgent(i) = 10.^((23-30)/10); % Case of maximum power transmission
    end
    
    for i=1:fbsNum
        for j=1:fbsNum
            if i ~= j
                if small_scale
                    P_interf = P_interf + pAgent(j)*(G(j,i)/L(j,i));
                else
                    P_interf = P_interf + pAgent(j)*(1/L(j,i));
                end
                
            end
        end
%         SINR(i) = (pAgent(i)*(G(i,i)/L(i,i)))/(MBS_P*(G(fbsNum+1,i)/L(fbsNum+1,i))+P_interf+sigma);
                if small_scale
                    SINR(i) = (pAgent(i)*(G(i,i)/L(i,i)))/(P_interf+sigma);
                else
                    SINR(i) = (pAgent(i)*(1/L(i,i)))/(P_interf+sigma);
                end        

        %SINR(i) = (pAgent(i)*(G(i,i)/L(i,i)))/(P_interf+sigma); % case when no MBS exists
        rx_inter(i)=P_interf;
    end
%     
    
%     for i=1:fbsNum  % In our case, UAVs are much apart, we can neglect their interferences.
%                 SINR(i) = (pAgent(i)*(G(i,i)/L(i,i)))/(MBS_P*(G(fbsNum+1,i)/L(fbsNum+1,i))+P_interf+sigma);
        %SINR(i) = (pAgent(i)*(G(i,i)/L(i,i)))/(P_interf+sigma); % case when no MBS exists
%     end
end