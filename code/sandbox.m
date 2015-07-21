clc
E1=1.23;
E12=3.21;
E3=2.13;
NM=1.42;

N=10;
ZE3=E3;

(NM^2/(2*E3^3))*(2*(E12+E3)*(expl(2,NM*E3)+expl(2,-NM*E3))-2*NM*E12*E3*sinh(NM*E3))

(NM^2/(2*E1^3))*(2*(E12+E1)*(expl(2,NM*E1)+expl(2,-NM*E1))-2*NM*E12*E1*sinh(NM*E1))

(NM^2/(E12^3))*(4*(E1+E12+E3)*sinh(0.5*NM*E12)^2-NM*(E1+E3)*E12*sinh(NM*E12))