

% E1=[1,1,1]'/(sqrt(3));
% E2=[1,1,1]'/(sqrt(3));

% E1=[1;0;0];
% E2=[1;0;0];


E1=[0;1;0];
E2=[0;1;0.000000001];


B1=acos(E1(3));
A1=atan2(E1(2),E1(1));

ROTZ=[cos(-A1) -sin(-A1) 0;sin(-A1) cos(-A1) 0;0 0 1];
ROTY=[cos(-B1) 0 sin(-B1);0 1 0;-sin(-B1) 0 cos(-B1)];

E2=ROTY*ROTZ*E2;

RHO=E2(3)
PHI=atan2(E2(2),E2(1))


%{
format long

L=1;
M=1;
theta=0.34748;
phi=0.38934748;

-0.5*sqrt(3/(2*pi()))*sin(theta)*exp(1i*phi)


RHO112=cos(theta);
PHI112=phi;
PLM112=legendre(L,RHO112);
%YLM112(L+1,M+1)=
exp(1i*M*PHI112)*sqrt(factorial(L-M)/factorial(L+M))*PLM112(M+1)
%}