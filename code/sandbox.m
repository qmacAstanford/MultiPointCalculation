function sandbox()

NM=1.4;
E1=0.5976;
E12=0.4357;
E3=0.7486;

NM^4*(-NM*(E1-2*E12+E3)+6)/12

NM*(2*(E1+E12+E3)*(expl(1,NM*E1)+expl(1,-NM*E1))...
    +NM*E1*(2*E1+NM*E1*E12+2*E12+E3*(exp(NM*E1)+1))*expl(1,-NM*E1))...
    /(2*E1^4)

NM*(2*(E1+E12+E3)*(expl(1,NM*E3)+expl(1,-NM*E3))...
    +NM*E3*(2*E3+NM*E3*E12+2*E12+E1*(exp(NM*E3)+1))*expl(1,-NM*E3))...
    /(2*E3^4)


NM^2*(2*(E1+E12+E3)*expl(2,NM*E12)-NM*(E1+E3)*E12*expl(1,NM*E12))/(2*E12^3)

% a=0.846546;
% b=0.75786;
% c=0.989753857;
% N=5;
% -c*(N*(b+a-5)+1-(N*(N+1)*((N^2-7*N+18)*(b+a-2)+4*(N-7)))/24)
% -c*b*(3*N-4*a+(N+1)*(18*a-3*N+N^2*a-7*N*a)/6-3)


% a=-0.342;
% b=0.1+0.4i;
% c=-0.2-0.3i;
% 
% abc=[a,b,c];
% [junk,Index]=max(real(abc));
% a=abc(Index);
% abc(Index)=[];
% [junk,Index]=min(abs(abc-a));
% b=abc(Index);
% abc(Index)=[];
% c=abc;
% 
% [a,b,c]




% E1= -78.3952-8.349*1i;
% E12=-26.5017-4.0378*1i;
% E1=-78+0.1i;
% E12=-26-2.1i;
% NM=10;
% 
% valeq=(exp(1).^E12).^((-1).*NM).*log(exp(1).^E1).^(-2).*(log(exp(1).^E1) ...
% +(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^(-2).*(((-1)+( ...
% exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*((-1)+(exp(1).^E1).^NM).* ...
% log(exp(1).^E12)).^2
% 
% 
% valeq=(E1^2*2*coshl(2,-E12*NM)+2*E1*E12*expl(1,E1*NM)*expl(1,-E12*NM)+...
%       E12^2*exp(-NM*E12)*expl(1,NM*E1)^2 )/((E1*E12*(E1-E12))^2)
%   
%   
  
end