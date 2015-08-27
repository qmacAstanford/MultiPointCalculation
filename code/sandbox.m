function sandbox()
NM=1.4;
E1=0.54;
E12=0.48;
E3=0.38;
addpath('Cases')


expl(1,NM*E1)*(expl(1,-NM*E1)-expl(1,-NM*E12))*...
              ( expl(2,NM*E3)*E12 - expl(2,NM*E12)*E3 )/...
              (E1*E12*E3*(E12-E3)*(E1-E12))

          
          
% % 0~E1
% valeq=(expl(1,NM*E3)*(expl(2,-NM*E12)+NM*E12*expl(1,-NM*E12))...
%     -expl(1,-NM*E3)*expl(2,NM*E12))/(E12^2*E3*(E12-E3))
% 
% % 0~E12
% valeq=-expl(1,NM*E3)*expl(1,-NM*E3)*expl(2,NM*E1)/(E1^2*E3^2)
% 
% % 0~E3
% valeq=NM*(E1*expl(2,NM*E12)+E1*expl(2,-NM*E12)...
%           +E12*expl(1,NM*E1)*expl(1,-NM*E12))...
%           /(E1*E12^2*(E12-E1))
%   
% % 0~E1, E12~E3
% valeq=-NM*expl(1,-NM*E12)*expl(2,NM*E12)/(E12^3)
% 
% % 0~E12, E1~E3
% valeq=(expl(2,NM*E1)+expl(2,-NM*E1))*expl(2,NM*E1)/(E1^4)
% 
% % 0~E3, E1~E12
% valeq=NM^2*expl(1,NM*E1)/(E1^2) - 4*NM*sinh(0.5*NM*E1)^2/(E1^3)
end