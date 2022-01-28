function d = nearest_FUE(fbs,FBS,fbsNum, i)
    d = inf;
    %d_hor_next= zeros(1,fbsNum-1);
    xAgent = fbs.X;
    yAgent = fbs.Y;
    hAgent = fbs.h;
    for j=1:fbsNum
        if i ~= j
            d_hor_next = sqrt((xAgent-FBS{j}.FUEX)^2+(yAgent-FBS{j}.FUEY)^2);
            d_next = sqrt((d_hor_next)^2+(hAgent)^2);  
            if d_next < d
                d = d_next;
            end
        end
    
    end
end