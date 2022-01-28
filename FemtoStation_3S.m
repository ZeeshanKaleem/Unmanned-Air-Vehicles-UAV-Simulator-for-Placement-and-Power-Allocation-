classdef FemtoStation_3S
   properties
      X
      Y
      h
      P
      dBS
      dMUE
      dFUE
      FUEX
      FUEY
      M  % distance with MUE
      B  % distance with BS
      eta_LOS=0.1; % Suburban
      eta_NLOS=21;% Suburban
%       eta_LOS=1.6; % dense urban
%       eta_NLOS=23;% dense urban
      a=4.88; % suburban scenario 
      b=0.43;% suburban scenario
      dM1 = 15; dM2 = 50; dM3 = 125; 
      %dB1 = 50; dB2 = 150; dB3 = 400;
      h1 = 100; h2 = 200; h3 = 300; h4 = 400; h5 = 500; h6 = 600; h7 = 700;
      state = zeros(1,2)
      powerProfile = []
      C_FUE
      d_near
      C_profile = []
      SINR_UAVUE
      SINR_profile = []
      R_FUE
      R_profile = []
      
   end
   methods
      function obj = FemtoStation_3S(xPos, yPos, zPos, BS, MUE, dFUE, rad_Ubs)
        obj.X = xPos;
        obj.Y = yPos;
        obj.h = zPos;
        obj.dBS = sqrt((xPos-BS.X)^2 + (yPos-BS.Y)^2);
        obj.dMUE = nearest_MUE(xPos, yPos, MUE);% sqrt((xPos-MUE.X)^2 + (yPos-MUE.Y)^2); %distance to nearest MUE
        obj.dFUE = dFUE;
%         obj.FUEX = xPos;
%         obj.FUEY = yPos+dFUE;
        obj.FUEX=(xPos-rad_Ubs) + ((xPos+rad_Ubs)-(xPos-rad_Ubs))*rand(1,1);
        obj.FUEY=(yPos-rad_Ubs) + ((yPos+rad_Ubs)-(yPos-rad_Ubs))*rand(1,1);
        
        
      end
      
      function obj = setPower(obj,power)
%           obj.P = 10^((power-30)/10);
            obj.P = power;
            obj.powerProfile = [obj.powerProfile power];
      end
      
      
      function obj = eSearchPower(obj,power)
%           obj.P = 10^((power-30)/10);
            obj.Pow = power;
            %obj.powerProfile = [obj.powerProfile power];
      end
      
      function obj = setCapacity(obj,c)
        obj.C_FUE = c;
        obj.C_profile = [obj.C_profile c];
      end
      
       function obj = reward(obj,r)
        obj.R_FUE = r;
        obj.R_profile = [obj.R_profile r];
      end
      
      function obj = SINR(obj,sinr)
        obj.SINR_UAVUE = 10*log10(sinr);
        obj.SINR_profile = [obj.SINR_profile 10*log10(sinr)];
      end
      
      function obj = min_neighbordistance(obj,near)
%           obj.P = 10^((power-30)/10);
            obj.d_near = near;
      end
      
      function obj = getDistanceStatus(obj, range, NUAVs,scenario)
         %here we are considering the different states in the x,y direction
         %just with different signs of the coordinates.
       if NUAVs==4  
         if(obj.X > 0 && obj.Y >0)
              obj.state(1) = 0;
              obj.state(2) = 0;
          elseif (obj.X < 0 && obj.Y >0)
              obj.state(1) = 1;
              obj.state(2) = 0;
          elseif (obj.X < 0 && obj.Y <0)
              obj.state(1) = 1;
              obj.state(2) = 1;
          else
              obj.state(1) = 0;
              obj.state(2) = 1;
         end
       elseif NUAVs==8
           if((0<obj.X && obj.X<= range) && (obj.Y>0 && obj.Y<2*range))
              obj.state(1) = 0;
              obj.state(2) = 0;
              obj.state(3) = 0;
          elseif ((range<obj.X && obj.X< 2*range) && (obj.Y>0 && obj.Y<2*range))
              obj.state(1) = 0;
              obj.state(2) = 0;
              obj.state(3) = 1;
          elseif ((obj.X > -range && obj.X<0) && (obj.Y>0 && obj.Y<2*range))
              obj.state(1) = 0;
              obj.state(2) = 1;
              obj.state(3) = 0;
          elseif ((-2*range<obj.X && obj.X< -range) && (obj.Y>0 && obj.Y<2*range))
              obj.state(1) = 1;
              obj.state(2) = 0;
              obj.state(3) = 0;
          elseif ((0<obj.X && obj.X<= range)  && (obj.Y<0 && obj.Y>-2*range))
              obj.state(1) = 1;
              obj.state(2) = 1;
              obj.state(3) = 0;
          elseif ((range<obj.X && obj.X< 2*range) && (obj.Y<0 && obj.Y>-2*range))
              obj.state(1) = 0;
              obj.state(2) = 1;
              obj.state(3) = 1;
          elseif ((obj.X > -range && obj.X<0) && (obj.Y<0 && obj.Y>-2*range))
              obj.state(1) = 1;
              obj.state(2) = 0;
              obj.state(3) = 1;
          else
              obj.state(1) = 1;
              obj.state(2) = 1;
              obj.state(3) = 1;
           end
         
           elseif NUAVs==16
              
               r=range;% 25
                r1=r+range; % 50
                r2 = r1+range; % 75
                r3 = r2+range; % 100
           
           if((0<obj.X && obj.X<= r) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 0;
              obj.state(3) = 0;
              obj.state(4) = 0;
          elseif ((r<obj.X && obj.X< r1) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 0;
              obj.state(3) = 0;
              obj.state(4) = 1;
          elseif ((r1<obj.X && obj.X< r2) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 0;
              obj.state(3) = 1;
              obj.state(4) = 0;
          elseif ((r2<obj.X && obj.X< r3) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 0;
              obj.state(3) = 1;
              obj.state(4) = 1;
          elseif ((obj.X > -r && obj.X<0) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 1;
              obj.state(3) = 0;
              obj.state(4) = 0;
          elseif ((-r1<obj.X && obj.X< -r) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 1;
              obj.state(3) = 0;
              obj.state(4) = 1;
          elseif ((-r2<obj.X && obj.X< -r1) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 1;
              obj.state(3) = 1;
              obj.state(4) = 0;
          elseif ((-r3<obj.X && obj.X< -r2) && (obj.Y>0 && obj.Y<4*r))
              obj.state(1) = 0;
              obj.state(2) = 1;
              obj.state(3) = 1;
              obj.state(4) = 1;
              
              
          elseif ((obj.X > -r && obj.X<0)  && (obj.Y<0 && obj.Y>-4*r))
              obj.state(1) = 1;
              obj.state(2) = 0;
              obj.state(3) = 0;
              obj.state(4) = 0;
          elseif ((-r1<obj.X && obj.X< -r) && (obj.Y<0 && obj.Y>-4*r))
              obj.state(1) = 1;
              obj.state(2) = 0;
              obj.state(3) = 0;
              obj.state(4) = 1;
          elseif ((-r2<obj.X && obj.X< -r1)  && (obj.Y<0 && obj.Y>-4*r))
              obj.state(1) = 1;
              obj.state(2) = 0;
              obj.state(3) = 1;
              obj.state(4) = 0;
          elseif ((-r3<obj.X && obj.X< -r2) && (obj.Y<0 && obj.Y>-4*r))
              obj.state(1) = 1;
              obj.state(2) = 0;
              obj.state(3) = 1;
              obj.state(4) = 1;
          elseif ((0<obj.X && obj.X<= r) && (obj.Y<0 && obj.Y>-4*r))
              obj.state(1) = 1;
              obj.state(2) = 1;
              obj.state(3) = 0;
              obj.state(4) = 0;
          elseif ((r<obj.X && obj.X< r1) && (obj.Y<0 && obj.Y>-4*r))
              obj.state(1) = 1;
              obj.state(2) = 1;
              obj.state(3) = 0;
              obj.state(4) = 1;
          elseif ((r1<obj.X && obj.X< r2) && (obj.Y<0 && obj.Y>-4*r))
              obj.state(1) = 1;
              obj.state(2) = 1;
              obj.state(3) = 1;
              obj.state(4) = 0;
          else
              obj.state(1) = 1;
              obj.state(2) = 1;
              obj.state(3) = 1;
              obj.state(4) = 1;
         end
      end
%           elseif(obj.dMUE <= obj.dM2 )
%               obj.state(1) = 1;
%           elseif(obj.dMUE <= obj.dM3 )
%               obj.state(1) = 2;
%           elseif
% %               obj.state(1) = 0;
%           elseif
%           else
          %end
          
%           if(obj.dBS <= obj.dB1 )
%               obj.state(2) = 0;
%           elseif(obj.dBS <= obj.dB2 )
%               obj.state(2) = 1;
%           elseif(obj.dBS <= obj.dB3 )
%               obj.state(2) = 2;
%           else
%               obj.state(2) = 0;
%           end
          
%          if(obj.h <= obj.h1 )
%              obj.state(3) = 0;
%          elseif(obj.h <= obj.h2 )
%              obj.state(3) = 1;
%          elseif(obj.h <= obj.h3)
%              obj.state(3) = 2;
%          elseif (obj.h <= obj.h4) 
%              obj.state(3) = 3;
%          elseif(obj.h <= obj.h5 )
%              obj.state(3) = 4;
%          elseif(obj.h <= obj.h6)
%              obj.state(3) = 5;
%          else
%              obj.state(3) = 6;
%          end
      end
     % function obj = getDistanceStatus(obj)
    %      if(obj.dMUE <= obj.dM1 )
    %          obj.state(2) = 0;
    %      elseif(obj.dMUE <= obj.dM2 )
    %          obj.state(2) = 1;
     %     elseif(obj.dMUE <= obj.dM3 )
    %          obj.state(2) = 2;
    %      else
    %          obj.state(2) = 3;
    %      end
          
   %       if(obj.dBS <= obj.dB1 )
   %           obj.state(3) = 0;
   %       elseif(obj.dBS <= obj.dB2 )
   %           obj.state(3) = 1;
   %       elseif(obj.dBS <= obj.dB3 )
  %            obj.state(3) = 2;
   %       else
   %           obj.state(3) = 3;
   %       end
     % end
   end
end