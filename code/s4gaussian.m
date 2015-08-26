function S4=s4gaussian(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d)
%S4 is a four point correlation function
%For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA

addpath('Cases')

% Reset Qs to column vectors if entered as rows
if isrow(Q1)==1
    Q1=transpose(Q1);
    Q2=transpose(Q2);
    Q3=transpose(Q3);
    Q4=transpose(Q4);
end

% Begin calculation of s4
MIN=1e-6;
if sum(power(Q1+Q2+Q3+Q4,2)) > MIN
    return
end

% Evaluate the quantities for s4 calculation
FB=1-FA;
F=[FA,FB];

%{
S4=zeros(2,2,2,2);

digits(5);
DF=[1,-1];

Q1MAG=sqrt(sum(power(Q1,2)));
Q2MAG=sqrt(sum(power(Q2,2)));
Q3MAG=sqrt(sum(power(Q3,2)));
Q4MAG=sqrt(sum(power(Q4,2)));

R1=-Q1MAG*Q1MAG/(2*d);
R2=-Q2MAG*Q2MAG/(2*d);
R3=-Q3MAG*Q3MAG/(2*d);
R4=-Q4MAG*Q4MAG/(2*d);

Q12MAG=sqrt(sum(power(Q1+Q2,2)));
Q13MAG=sqrt(sum(power(Q1+Q3,2)));
Q14MAG=sqrt(sum(power(Q1+Q4,2)));
Q23MAG=sqrt(sum(power(Q2+Q3,2)));
Q24MAG=sqrt(sum(power(Q2+Q4,2)));
Q34MAG=sqrt(sum(power(Q3+Q4,2)));

R12=-Q12MAG*Q12MAG/(2*d);
R13=-Q13MAG*Q13MAG/(2*d);
R14=-Q14MAG*Q14MAG/(2*d);
R23=-Q23MAG*Q23MAG/(2*d);
R24=-Q24MAG*Q24MAG/(2*d);
R34=-Q34MAG*Q34MAG/(2*d);

Z1=exp(R1*NM);
Z2=exp(R2*NM);
Z3=exp(R3*NM);
Z4=exp(R4*NM);

Z12=exp(R12*NM);
Z13=exp(R13*NM);
Z14=exp(R14*NM);
Z23=exp(R23*NM);
Z24=exp(R24*NM);
Z34=exp(R34*NM);

Z1L=exp(R1*NM)*LAM;
Z2L=exp(R2*NM)*LAM;
Z3L=exp(R3*NM)*LAM;
Z4L=exp(R4*NM)*LAM;

Z12L=exp(R12*NM)*LAM;
Z14L=exp(R14*NM)*LAM;
Z13L=exp(R13*NM)*LAM;
Z23L=exp(R23*NM)*LAM;
Z24L=exp(R24*NM)*LAM;
Z34L=exp(R34*NM)*LAM;

% PERMUTATIONS
% 1 12 4 112 443
% 1 12 3 112 334
% 1 13 4 113 442
% 1 13 2 113 224
% 1 14 3 114 332
% 1 14 2 114 223
% 2 23 4 223 441
% 2 23 1 223 114
% 2 24 3 224 331
% 2 24 1 224 113
% 3 34 2 334 221
% 3 34 1 334 112
CasesIncluded=[0,0,0,0,0,0,1,0];
%              1 2 3 4 5 6 7 8

if CasesIncluded(1)
% Case 1: J1==J2==J3==J4
    S4=case1(S4,N,NM,R1,R12,R4,F,MIN);
    S4=case1(S4,N,NM,R1,R12,R3,F,MIN);
    S4=case1(S4,N,NM,R1,R13,R4,F,MIN);
    S4=case1(S4,N,NM,R1,R13,R2,F,MIN);
    S4=case1(S4,N,NM,R1,R14,R3,F,MIN);
    S4=case1(S4,N,NM,R1,R14,R2,F,MIN);
    S4=case1(S4,N,NM,R2,R23,R4,F,MIN);
    S4=case1(S4,N,NM,R2,R23,R1,F,MIN);
    S4=case1(S4,N,NM,R2,R24,R3,F,MIN);
    S4=case1(S4,N,NM,R2,R24,R1,F,MIN);
    S4=case1(S4,N,NM,R3,R34,R2,F,MIN);
    S4=case1(S4,N,NM,R3,R34,R1,F,MIN);
end
if CasesIncluded(2)
% Case 2: J1 <J2==J3==J4 (4 subcases)
    % J1 <J2==J3==J4 (six terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,2,2,2)=1;
    %total of six terms
    S4=case2(S4,N,NM,R1,R12,R4,Z1,Z1L,F,DF,SDEL,1,2,MIN);
    S4=case2(S4,N,NM,R1,R12,R3,Z1,Z1L,F,DF,SDEL,1,2,MIN);
    S4=case2(S4,N,NM,R1,R13,R4,Z1,Z1L,F,DF,SDEL,1,3,MIN);
    S4=case2(S4,N,NM,R1,R13,R2,Z1,Z1L,F,DF,SDEL,1,3,MIN);
    S4=case2(S4,N,NM,R1,R14,R3,Z1,Z1L,F,DF,SDEL,1,4,MIN);
    S4=case2(S4,N,NM,R1,R14,R2,Z1,Z1L,F,DF,SDEL,1,4,MIN);

    % J2 <J1==J3==J4
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,2,1,1)=1;
    SDEL(2,1,2,2)=1;
    S4=case2(S4,N,NM,R2,R12,R3,Z2,Z2L,F,DF,SDEL,2,1,MIN);
    S4=case2(S4,N,NM,R2,R12,R4,Z2,Z2L,F,DF,SDEL,2,1,MIN);
    S4=case2(S4,N,NM,R2,R23,R1,Z2,Z2L,F,DF,SDEL,2,3,MIN);
    S4=case2(S4,N,NM,R2,R23,R4,Z2,Z2L,F,DF,SDEL,2,3,MIN);
    S4=case2(S4,N,NM,R2,R24,R1,Z2,Z2L,F,DF,SDEL,2,4,MIN);
    S4=case2(S4,N,NM,R2,R24,R3,Z2,Z2L,F,DF,SDEL,2,4,MIN);

    % J3 <J1==J2==J4
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,2,1)=1;
    SDEL(2,2,1,2)=1;
    S4=case2(S4,N,NM,R3,R13,R2,Z2,Z2L,F,DF,SDEL,3,1,MIN);
    S4=case2(S4,N,NM,R3,R13,R4,Z2,Z2L,F,DF,SDEL,3,1,MIN);
    S4=case2(S4,N,NM,R3,R23,R1,Z2,Z2L,F,DF,SDEL,3,2,MIN);
    S4=case2(S4,N,NM,R3,R23,R4,Z2,Z2L,F,DF,SDEL,3,2,MIN);
    S4=case2(S4,N,NM,R3,R34,R1,Z2,Z2L,F,DF,SDEL,3,4,MIN);
    S4=case2(S4,N,NM,R3,R34,R2,Z2,Z2L,F,DF,SDEL,3,4,MIN);

    % J4 <J1==J2==J3
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,1,2)=1;
    SDEL(2,2,2,1)=1;
    S4=case2(S4,N,NM,R4,R14,R2,Z2,Z2L,F,DF,SDEL,4,1,MIN);
    S4=case2(S4,N,NM,R4,R14,R3,Z2,Z2L,F,DF,SDEL,4,1,MIN);
    S4=case2(S4,N,NM,R4,R24,R1,Z2,Z2L,F,DF,SDEL,4,2,MIN);
    S4=case2(S4,N,NM,R4,R24,R3,Z2,Z2L,F,DF,SDEL,4,2,MIN);
    S4=case2(S4,N,NM,R4,R34,R1,Z2,Z2L,F,DF,SDEL,4,3,MIN);
    S4=case2(S4,N,NM,R4,R34,R2,Z2,Z2L,F,DF,SDEL,4,3,MIN);
end
if CasesIncluded(3)
% Case 3: J1==J2 <J3==J4 (3 subcases)
    % J1=J2<J3=J4 // J2=J1<J3=J4 // J1=J2<J4=J3 // J2=J1<J4=J3
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,2,2)=1;
    SDEL(2,2,1,1)=1;
    S4=case3(S4,N,NM,R1,R12,R4,Z12,Z12L,F,DF,SDEL,2,3,MIN);
    S4=case3(S4,N,NM,R2,R12,R4,Z12,Z12L,F,DF,SDEL,1,3,MIN);
    S4=case3(S4,N,NM,R1,R12,R3,Z12,Z12L,F,DF,SDEL,2,4,MIN);
    S4=case3(S4,N,NM,R2,R12,R3,Z12,Z12L,F,DF,SDEL,1,4,MIN);
    % J3=J4<J1=J2 // J4=J3<J1=J2 // J3=J4<J2=J1 // J4=J3<J2=J1
    S4=case3(S4,N,NM,R3,R34,R2,Z34,Z34L,F,DF,SDEL,4,1,MIN);
    S4=case3(S4,N,NM,R4,R34,R2,Z34,Z34L,F,DF,SDEL,3,1,MIN);
    S4=case3(S4,N,NM,R3,R34,R1,Z34,Z34L,F,DF,SDEL,4,2,MIN);
    S4=case3(S4,N,NM,R4,R34,R1,Z34,Z34L,F,DF,SDEL,3,2,MIN);

    % J1=J3<J2=J4 // J3=J1<J2=J4 // J1=J3<J4=J2 // J3=J1<J4=J2
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,2,1,2)=1;
    SDEL(2,1,2,1)=1;
    S4=case3(S4,N,NM,R1,R13,R4,Z13,Z13L,F,DF,SDEL,3,2,MIN);
    S4=case3(S4,N,NM,R3,R13,R4,Z13,Z13L,F,DF,SDEL,1,2,MIN);
    S4=case3(S4,N,NM,R1,R13,R2,Z13,Z13L,F,DF,SDEL,3,4,MIN);
    S4=case3(S4,N,NM,R3,R13,R2,Z13,Z13L,F,DF,SDEL,1,4,MIN);
    % J2=J4<J1=J3 // J2=J4<J3=J1 // J4=J2<J1=J3 // J4=J2<J3=J1
    S4=case3(S4,N,NM,R2,R24,R3,Z24,Z24L,F,DF,SDEL,4,1,MIN);
    S4=case3(S4,N,NM,R2,R24,R1,Z24,Z24L,F,DF,SDEL,4,3,MIN);
    S4=case3(S4,N,NM,R4,R24,R3,Z24,Z24L,F,DF,SDEL,2,1,MIN);
    S4=case3(S4,N,NM,R4,R24,R1,Z24,Z24L,F,DF,SDEL,2,3,MIN);

    % J1=J4<J2=J3 // J4=J1<J2=J3 // J1=J4<J3=J2 // J4=J1<J3=J2
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,2,2,1)=1;
    SDEL(2,1,1,2)=1;
    % J2=J3<J1=J4 // J3=J2<J1=J4 // J2=J3<J4=J1 // J3=J2<J4=J1
    S4=case3(S4,N,NM,R2,R23,R4,Z23,Z23L,F,DF,SDEL,3,1,MIN);
    S4=case3(S4,N,NM,R3,R23,R4,Z23,Z23L,F,DF,SDEL,2,1,MIN);
    S4=case3(S4,N,NM,R2,R23,R1,Z23,Z23L,F,DF,SDEL,3,4,MIN);
    S4=case3(S4,N,NM,R3,R23,R1,Z23,Z23L,F,DF,SDEL,2,4,MIN);
    % J1=J4<J2=J3 // J1=J4<J3=J2 // J4=J1<J2=J3 // J4=J1<J3=J2
    S4=case3(S4,N,NM,R1,R14,R3,Z14,Z14L,F,DF,SDEL,4,2,MIN);
    S4=case3(S4,N,NM,R1,R14,R2,Z14,Z14L,F,DF,SDEL,4,3,MIN);
    S4=case3(S4,N,NM,R4,R14,R3,Z14,Z14L,F,DF,SDEL,1,2,MIN);
    S4=case3(S4,N,NM,R4,R14,R2,Z14,Z14L,F,DF,SDEL,1,3,MIN);
end
if CasesIncluded(4)
% Case 4: J1==J2==J3 <J4 (4 subcases)
    % J1==J2==J3 <J4 (six terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,1,2)=1;
    SDEL(2,2,2,1)=1;
    S4=case4(S4,N,NM,R1,R12,R4,Z4,Z4L,F,DF,SDEL,3,4,MIN);
    S4=case4(S4,N,NM,R1,R13,R4,Z4,Z4L,F,DF,SDEL,2,4,MIN);
    S4=case4(S4,N,NM,R2,R12,R4,Z4,Z4L,F,DF,SDEL,3,4,MIN);
    S4=case4(S4,N,NM,R2,R23,R4,Z4,Z4L,F,DF,SDEL,1,4,MIN);
    S4=case4(S4,N,NM,R3,R13,R4,Z4,Z4L,F,DF,SDEL,2,4,MIN);
    S4=case4(S4,N,NM,R3,R23,R4,Z4,Z4L,F,DF,SDEL,1,4,MIN);    

    % J1==J2==J4 <J3
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,2,1)=1;
    SDEL(2,2,1,2)=1;
    S4=case4(S4,N,NM,R1,R12,R3,Z3,Z3L,F,DF,SDEL,4,3,MIN);
    S4=case4(S4,N,NM,R1,R14,R3,Z3,Z3L,F,DF,SDEL,2,3,MIN);
    S4=case4(S4,N,NM,R2,R12,R3,Z3,Z3L,F,DF,SDEL,4,3,MIN);
    S4=case4(S4,N,NM,R2,R24,R3,Z3,Z3L,F,DF,SDEL,1,3,MIN);
    S4=case4(S4,N,NM,R4,R14,R3,Z3,Z3L,F,DF,SDEL,2,3,MIN);
    S4=case4(S4,N,NM,R4,R24,R3,Z3,Z3L,F,DF,SDEL,1,3,MIN);

    % J1==J3==J4 <J2
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,2,1,1)=1;
    SDEL(2,1,2,2)=1;
    S4=case4(S4,N,NM,R1,R13,R2,Z2,Z2L,F,DF,SDEL,4,2,MIN);
    S4=case4(S4,N,NM,R1,R14,R2,Z2,Z2L,F,DF,SDEL,3,2,MIN);
    S4=case4(S4,N,NM,R3,R13,R2,Z2,Z2L,F,DF,SDEL,4,2,MIN);
    S4=case4(S4,N,NM,R3,R34,R2,Z2,Z2L,F,DF,SDEL,1,2,MIN);
    S4=case4(S4,N,NM,R4,R14,R2,Z2,Z2L,F,DF,SDEL,3,2,MIN);
    S4=case4(S4,N,NM,R4,R34,R2,Z2,Z2L,F,DF,SDEL,1,2,MIN);

    % J2==J3==J4 <J1
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,2,2,2)=1;
    S4=case4(S4,N,NM,R2,R23,R1,Z1,Z1L,F,DF,SDEL,4,1,MIN);
    S4=case4(S4,N,NM,R2,R24,R1,Z1,Z1L,F,DF,SDEL,3,1,MIN);
    S4=case4(S4,N,NM,R3,R23,R1,Z1,Z1L,F,DF,SDEL,4,1,MIN);
    S4=case4(S4,N,NM,R3,R34,R1,Z1,Z1L,F,DF,SDEL,2,1,MIN);
    S4=case4(S4,N,NM,R4,R24,R1,Z1,Z1L,F,DF,SDEL,3,1,MIN);
    S4=case4(S4,N,NM,R4,R34,R1,Z1,Z1L,F,DF,SDEL,2,1,MIN);
end
if CasesIncluded(5)
% Case 5: J1 <J2 <J3==J4 (six subcases)
    % J1 <J2 <J3==J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,2,1,1)=1;
    SDEL(2,2,1,1)=1;
    SDEL(1,2,2,2)=1;
    SDEL(2,1,2,2)=1;
    SDEL(1,1,2,2)=1;
    S4=case5(S4,N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,F,DF,SDEL,1,2,3,MIN);
    S4=case5(S4,N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,F,DF,SDEL,1,2,4,MIN);
    % J2 <J1 <J3==J4
    S4=case5(S4,N,NM,R2,R12,R4,Z1,Z2L,Z12,Z12L,F,DF,SDEL,2,1,3,MIN);
    S4=case5(S4,N,NM,R2,R12,R3,Z1,Z2L,Z12,Z12L,F,DF,SDEL,2,1,4,MIN);

    % J1 <J3 <J2==J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,1,2,1)=1;
    SDEL(2,1,2,1)=1;
    SDEL(1,2,2,2)=1;
    SDEL(2,2,1,2)=1;
    SDEL(1,2,1,2)=1;
    S4=case5(S4,N,NM,R1,R13,R4,Z1,Z1L,Z13,Z13L,F,DF,SDEL,1,3,2,MIN);
    S4=case5(S4,N,NM,R1,R13,R2,Z1,Z1L,Z13,Z13L,F,DF,SDEL,1,3,4,MIN);
    % J3 <J1 <J2==J4
    S4=case5(S4,N,NM,R3,R13,R4,Z3,Z3L,Z13,Z13L,F,DF,SDEL,3,1,2,MIN);
    S4=case5(S4,N,NM,R3,R13,R2,Z3,Z3L,Z13,Z13L,F,DF,SDEL,3,1,4,MIN);

    % J2 <J3 <J1==J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,2,1)=1;
    SDEL(1,2,1,1)=1;
    SDEL(1,2,2,1)=1;
    SDEL(2,1,2,2)=1;
    SDEL(2,2,1,2)=1;
    SDEL(2,1,1,2)=1;
    S4=case5(S4,N,NM,R2,R23,R4,Z2,Z2L,Z23,Z23L,F,DF,SDEL,2,3,1,MIN);
    S4=case5(S4,N,NM,R2,R23,R1,Z2,Z2L,Z23,Z23L,F,DF,SDEL,2,3,4,MIN);
    % J3 <J2 <J1==J4
    S4=case5(S4,N,NM,R3,R23,R4,Z3,Z3L,Z23,Z23L,F,DF,SDEL,3,2,1,MIN);
    S4=case5(S4,N,NM,R3,R23,R1,Z3,Z3L,Z23,Z23L,F,DF,SDEL,3,2,4,MIN);

    % J1 <J4 <J2==J3 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,1,1,2)=1;
    SDEL(2,1,1,2)=1;
    SDEL(1,2,2,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(1,2,2,1)=1;
    S4=case5(S4,N,NM,R1,R14,R3,Z1,Z1L,Z14,Z14L,F,DF,SDEL,1,4,2,MIN);
    S4=case5(S4,N,NM,R1,R14,R2,Z1,Z1L,Z14,Z14L,F,DF,SDEL,1,4,3,MIN);
    % J4 <J1 <J2==J3
    S4=case5(S4,N,NM,R4,R14,R3,Z4,Z4L,Z14,Z14L,F,DF,SDEL,4,1,2,MIN);
    S4=case5(S4,N,NM,R4,R14,R2,Z4,Z4L,Z14,Z14L,F,DF,SDEL,4,1,3,MIN);

    % J2 <J4 <J1==J3 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,2,1,1)=1;
    SDEL(1,1,1,2)=1;
    SDEL(1,2,1,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(2,1,2,2)=1;
    SDEL(2,1,2,1)=1;
    S4=case5(S4,N,NM,R2,R24,R3,Z2,Z2L,Z24,Z24L,F,DF,SDEL,2,4,1,MIN);
    S4=case5(S4,N,NM,R2,R24,R1,Z2,Z2L,Z24,Z24L,F,DF,SDEL,2,4,3,MIN);
    % J4 <J2 <J1==J3
    S4=case5(S4,N,NM,R4,R24,R3,Z4,Z4L,Z24,Z24L,F,DF,SDEL,4,2,1,MIN);
    S4=case5(S4,N,NM,R4,R24,R1,Z4,Z4L,Z24,Z24L,F,DF,SDEL,4,2,3,MIN);

    % J3 <J4 <J1==J2 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,2,1)=1;
    SDEL(1,1,1,2)=1;
    SDEL(1,1,2,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(2,2,1,2)=1;
    SDEL(2,2,1,1)=1;
    S4=case5(S4,N,NM,R3,R34,R2,Z3,Z3L,Z34,Z34L,F,DF,SDEL,3,4,1,MIN);
    S4=case5(S4,N,NM,R3,R34,R1,Z3,Z3L,Z34,Z34L,F,DF,SDEL,3,4,2,MIN);
    % J4 <J3 <J1==J2
    S4=case5(S4,N,NM,R4,R34,R2,Z4,Z4L,Z34,Z34L,F,DF,SDEL,4,3,1,MIN);
    S4=case5(S4,N,NM,R4,R34,R1,Z4,Z4L,Z34,Z34L,F,DF,SDEL,4,3,2,MIN);
end
if CasesIncluded(6)
% Case 6: J1==J2 <J3 <J4 (six subterms)
    % J1==J2 <J3 <J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,1,2)=1;
    SDEL(1,1,2,1)=1;
    SDEL(1,1,2,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(2,2,1,2)=1;
    SDEL(2,2,1,1)=1;
    S4=case6(S4,N,NM,R1,R12,R4,Z12,Z12L,Z4,Z4L,F,DF,SDEL,2,3,4,MIN);
    S4=case6(S4,N,NM,R2,R12,R4,Z12,Z12L,Z4,Z4L,F,DF,SDEL,1,3,4,MIN);        
    % J1==J2 <J4 <J3
    S4=case6(S4,N,NM,R1,R12,R3,Z12,Z12L,Z3,Z3L,F,DF,SDEL,2,4,3,MIN);
    S4=case6(S4,N,NM,R2,R12,R3,Z12,Z12L,Z3,Z3L,F,DF,SDEL,1,4,3,MIN);

    % J1==J3 <J2 <J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,1,2)=1;
    SDEL(1,2,1,1)=1;
    SDEL(1,2,1,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(2,1,2,2)=1;
    SDEL(2,1,2,1)=1;
    S4=case6(S4,N,NM,R1,R13,R4,Z13,Z13L,Z4,Z4L,F,DF,SDEL,3,2,4,MIN);
    S4=case6(S4,N,NM,R3,R13,R4,Z13,Z13L,Z4,Z4L,F,DF,SDEL,1,2,4,MIN);
    % J1==J3 <J4 <J2
    S4=case6(S4,N,NM,R1,R13,R2,Z13,Z13L,Z2,Z2L,F,DF,SDEL,3,4,2,MIN);
    S4=case6(S4,N,NM,R3,R13,R2,Z13,Z13L,Z2,Z2L,F,DF,SDEL,1,4,2,MIN);

    % J1==J4 <J2 <J3 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,2,1)=1;
    SDEL(1,2,1,1)=1;
    SDEL(1,2,2,1)=1;
    SDEL(2,2,1,2)=1;
    SDEL(2,1,2,2)=1;
    SDEL(2,1,1,2)=1;
    S4=case6(S4,N,NM,R1,R14,R3,Z14,Z14L,Z3,Z3L,F,DF,SDEL,4,2,3,MIN);
    S4=case6(S4,N,NM,R4,R14,R3,Z14,Z14L,Z3,Z3L,F,DF,SDEL,1,2,3,MIN);
    % J1==J4 <J3 <J2
    S4=case6(S4,N,NM,R1,R14,R2,Z14,Z14L,Z2,Z2L,F,DF,SDEL,4,3,2,MIN);
    S4=case6(S4,N,NM,R4,R14,R2,Z14,Z14L,Z2,Z2L,F,DF,SDEL,1,3,2,MIN);

    % J2==J3 <J1 <J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,1,1,2)=1;
    SDEL(2,1,1,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(1,2,2,2)=1;
    SDEL(1,2,2,1)=1;
    S4=case6(S4,N,NM,R2,R23,R4,Z23,Z23L,Z4,Z4L,F,DF,SDEL,3,1,4,MIN);
    S4=case6(S4,N,NM,R3,R23,R4,Z23,Z23L,Z4,Z4L,F,DF,SDEL,2,1,4,MIN);
    % J2==J3 <J4 <J1
    S4=case6(S4,N,NM,R2,R23,R1,Z23,Z23L,Z1,Z1L,F,DF,SDEL,3,4,1,MIN);
    S4=case6(S4,N,NM,R3,R23,R1,Z23,Z23L,Z1,Z1L,F,DF,SDEL,2,4,1,MIN);

    % J2==J4 <J1 <J3 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,1,2,1)=1;
    SDEL(2,1,2,1)=1;
    SDEL(1,2,2,2)=1;
    SDEL(2,2,1,2)=1;
    SDEL(1,2,1,2)=1;
    S4=case6(S4,N,NM,R2,R24,R3,Z24,Z24L,Z3,Z3L,F,DF,SDEL,2,1,3,MIN);
    S4=case6(S4,N,NM,R4,R24,R3,Z24,Z24L,Z3,Z3L,F,DF,SDEL,4,1,3,MIN);
    % J2==J4 <J3 <J1
    S4=case6(S4,N,NM,R2,R24,R1,Z24,Z24L,Z1,Z1L,F,DF,SDEL,4,3,1,MIN);
    S4=case6(S4,N,NM,R4,R24,R1,Z24,Z24L,Z1,Z1L,F,DF,SDEL,2,3,1,MIN);

    % J3==J4 <J1 <J2 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,2,1,1)=1;
    SDEL(2,2,1,1)=1;
    SDEL(1,2,2,2)=1;
    SDEL(2,1,2,2)=1;
    SDEL(1,1,2,2)=1;
    S4=case6(S4,N,NM,R3,R34,R2,Z34,Z34L,Z2,Z2L,F,DF,SDEL,4,1,2,MIN);
    S4=case6(S4,N,NM,R4,R34,R2,Z34,Z34L,Z2,Z2L,F,DF,SDEL,3,1,2,MIN);        
    % J3==J4 <J2 <J1
    S4=case6(S4,N,NM,R3,R34,R1,Z34,Z34L,Z1,Z1L,F,DF,SDEL,4,2,1,MIN);
    S4=case6(S4,N,NM,R4,R34,R1,Z34,Z34L,Z1,Z1L,F,DF,SDEL,3,2,1,MIN);
end
if CasesIncluded(7)
% Case 7: J1 <J2==J3 <J4 (six subcases)
    % J1 <J2==J3 <J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,1,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(2,1,1,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(1,2,2,2)=1;
    SDEL(1,2,2,1)=1;
    disp('old')
    S4=case7(S4,N,NM,R1,R12,R4,Z1,Z1L,Z4,Z4L,F,DF,SDEL,1,2,4,MIN);
    S4=case7(S4,N,NM,R1,R13,R4,Z1,Z1L,Z4,Z4L,F,DF,SDEL,1,3,4,MIN);
    % J4 <J2==J3 <J1
    S4=case7(S4,N,NM,R4,R12,R1,Z4,Z4L,Z1,Z1L,F,DF,SDEL,4,2,1,MIN);
    S4=case7(S4,N,NM,R4,R13,R1,Z4,Z4L,Z1,Z1L,F,DF,SDEL,4,3,1,MIN);

    % J1 <J2==J4 <J3 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(2,1,1,1)=1;
    SDEL(1,1,2,1)=1;
    SDEL(2,1,2,1)=1;
    SDEL(1,2,2,2)=1;
    SDEL(2,2,1,2)=1;
    SDEL(1,2,1,2)=1;
    S4=case7(S4,N,NM,R1,R12,R3,Z1,Z1L,Z3,Z3L,F,DF,SDEL,1,2,3,MIN);
    S4=case7(S4,N,NM,R1,R14,R3,Z1,Z1L,Z3,Z3L,F,DF,SDEL,1,4,3,MIN);
    % J3 <J2==J4 <J1
    S4=case7(S4,N,NM,R3,R23,R1,Z3,Z3L,Z1,Z1L,F,DF,SDEL,3,2,1,MIN);
    S4=case7(S4,N,NM,R3,R34,R1,Z3,Z3L,Z1,Z1L,F,DF,SDEL,3,4,1,MIN);

    % J3 <J2==J1 <J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,1,2)=1;
    SDEL(1,1,2,1)=1;
    SDEL(1,1,2,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(2,2,1,2)=1;
    SDEL(2,2,1,1)=1;
    S4=case7(S4,N,NM,R3,R23,R4,Z3,Z3L,Z4,Z4L,F,DF,SDEL,3,2,4,MIN);
    S4=case7(S4,N,NM,R3,R13,R4,Z3,Z3L,Z4,Z4L,F,DF,SDEL,3,1,4,MIN);
    % J4 <J2==J1 <J3
    S4=case7(S4,N,NM,R4,R24,R3,Z4,Z4L,Z3,Z3L,F,DF,SDEL,4,2,3,MIN);
    S4=case7(S4,N,NM,R4,R14,R3,Z4,Z4L,Z3,Z3L,F,DF,SDEL,4,1,3,MIN);

S4=S4*0;    
    % J2 <J1==J3 <J4 (two terms)
    SDEL=zeros(2,2,2,2);
    SDEL(1,1,1,1)=1;
    SDEL(2,2,2,2)=1;
    SDEL(1,1,1,2)=1;
    SDEL(1,2,1,1)=1;
    SDEL(1,2,1,2)=1;
    SDEL(2,2,2,1)=1;
    SDEL(2,1,2,2)=1;
    SDEL(2,1,2,1)=1;
    %S4=case7(S4,N,NM,R2,R12,R4,Z2,Z2L,Z4,Z4L,F,DF,SDEL,2,1,4,MIN);
    S4=case7(S4,N,NM,R2,R13,R4,Z2,Z2L,Z4,Z4L,F,DF,SDEL,2,3,4,MIN);
    % J4 <J1==J3 <J2
    %S4=case7(S4,N,NM,R4,R12,R2,Z4,Z4L,Z2,Z2L,F,DF,SDEL,4,1,2,MIN);
    %S4=case7(S4,N,NM,R4,R13,R2,Z4,Z4L,Z2,Z2L,F,DF,SDEL,4,3,2,MIN);

    
%     % J2 <J1==J4 <J3 (two terms)
%     SDEL=zeros(2,2,2,2);
%     SDEL(1,1,1,1)=1;
%     SDEL(2,2,2,2)=1;
%     SDEL(1,1,2,1)=1;
%     SDEL(1,2,1,1)=1;
%     SDEL(1,2,2,1)=1;
%     SDEL(2,1,2,2)=1;
%     SDEL(2,2,1,2)=1;
%     SDEL(2,1,1,2)=1;
%     S4=case7(S4,N,NM,R2,R12,R3,Z2,Z2L,Z3,Z3L,F,DF,SDEL,2,1,3,MIN);
%     S4=case7(S4,N,NM,R2,R24,R3,Z2,Z2L,Z3,Z3L,F,DF,SDEL,2,4,3,MIN);
%     % J3 <J1==J4 <J2
%     S4=case7(S4,N,NM,R3,R13,R2,Z3,Z3L,Z2,Z2L,F,DF,SDEL,3,1,2,MIN);
%     S4=case7(S4,N,NM,R3,R34,R2,Z3,Z3L,Z2,Z2L,F,DF,SDEL,3,4,2,MIN);

    
    
%     % J1 <J3==J4 <J2 (two terms)
%     SDEL=zeros(2,2,2,2);
%     SDEL(1,1,1,1)=1;
%     SDEL(2,2,2,2)=1;
%     SDEL(2,1,1,1)=1;
%     SDEL(1,2,1,1)=1;
%     SDEL(2,2,1,1)=1;
%     SDEL(1,2,2,2)=1;  
%     SDEL(2,1,2,2)=1;   
%     SDEL(1,1,2,2)=1;
%     S4=case7(S4,N,NM,R1,R13,R2,Z1,Z1L,Z2,Z2L,F,DF,SDEL,1,3,2,MIN);
%     S4=case7(S4,N,NM,R1,R14,R2,Z1,Z1L,Z2,Z2L,F,DF,SDEL,1,4,2,MIN);
%     % J2 <J3==J4 <J1
%     S4=case7(S4,N,NM,R2,R23,R1,Z2,Z2L,Z1,Z1L,F,DF,SDEL,2,3,1,MIN);
%     S4=case7(S4,N,NM,R2,R24,R1,Z2,Z2L,Z1,Z1L,F,DF,SDEL,2,4,1,MIN);
    
    disp('old:')
    disp(S4)
    S4=S4*0;
    
    orders = perms(1:4);
    for orderNum=1:24
        SDEL=ones(2,2,2,2);
        order=orders(orderNum,:);
        Q = [Q1,Q2,Q3,Q4];
        Qnew=[Q(:,order(1)),Q(:,order(2)),Q(:,order(3)),Q(:,order(4))];
        
        R1=-dot(Qnew(:,1),Qnew(:,1))/(2*d);
        R12=-dot(Qnew(:,1)+Qnew(:,2),Qnew(:,1)+Qnew(:,2))/(2*d);
        R4=-dot(Qnew(:,4),Qnew(:,4))/(2*d);
        Z1=exp(R1*NM);
        Z12=exp(R12*NM);
        Z4=exp(R4*NM);
        
        if isequal(order,[2,3,1,4]) % || isequal(order,[2,1,3,4]) || isequal(order,[4,1,3,2]) || isequal(order,[4,3,1,2])
            %disp('new')
            S4=case7v2(S4,N,NM,R1,R12,R4,Z1,Z1*LAM,Z4,Z4*LAM,F,NaN,NaN,order,NaN);
        end
        
        %S4=case7v2(S4,N,NM,R1,R12,R3,Z1,Z1*LAM,Z4,Z4*LAM,F,NaN,NaN,order,NaN);
     

    end    
    disp('new')
    disp(S4)
    
    
end

if CasesIncluded(8)
% Case 8: J1 <J2 <J3 <J4


    SDEL=ones(2,2,2,2);
    S4=case8(S4,N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z4,Z4L,F,DF,SDEL,1,2,3,4,MIN);
    S4=case8(S4,N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,DF,SDEL,1,2,4,3,MIN);
    S4=case8(S4,N,NM,R1,R13,R4,Z1,Z1L,Z13,Z13L,Z4,Z4L,F,DF,SDEL,1,3,2,4,MIN);
    S4=case8(S4,N,NM,R1,R13,R2,Z1,Z1L,Z13,Z13L,Z2,Z2L,F,DF,SDEL,1,3,4,2,MIN);
    S4=case8(S4,N,NM,R1,R14,R3,Z1,Z1L,Z14,Z14L,Z3,Z3L,F,DF,SDEL,1,4,2,3,MIN);
    S4=case8(S4,N,NM,R1,R14,R2,Z1,Z1L,Z14,Z14L,Z2,Z2L,F,DF,SDEL,1,4,3,2,MIN);
    S4=case8(S4,N,NM,R2,R23,R4,Z2,Z2L,Z23,Z23L,Z4,Z4L,F,DF,SDEL,2,3,1,4,MIN);
    S4=case8(S4,N,NM,R2,R23,R1,Z2,Z2L,Z23,Z23L,Z1,Z1L,F,DF,SDEL,2,3,4,1,MIN);
    S4=case8(S4,N,NM,R2,R24,R3,Z2,Z2L,Z24,Z24L,Z3,Z3L,F,DF,SDEL,2,4,1,3,MIN);
    S4=case8(S4,N,NM,R2,R24,R1,Z2,Z2L,Z24,Z24L,Z1,Z1L,F,DF,SDEL,2,4,3,1,MIN);
    S4=case8(S4,N,NM,R3,R34,R2,Z3,Z3L,Z34,Z34L,Z2,Z2L,F,DF,SDEL,3,4,1,2,MIN);
    S4=case8(S4,N,NM,R3,R34,R1,Z3,Z3L,Z34,Z34L,Z1,Z1L,F,DF,SDEL,3,4,2,1,MIN);

    S4=case8(S4,N,NM,R2,R12,R4,Z2,Z2L,Z12,Z12L,Z4,Z4L,F,DF,SDEL,2,1,3,4,MIN);
    S4=case8(S4,N,NM,R2,R12,R3,Z2,Z2L,Z12,Z12L,Z3,Z3L,F,DF,SDEL,2,1,4,3,MIN);
    S4=case8(S4,N,NM,R3,R13,R4,Z3,Z3L,Z13,Z13L,Z4,Z4L,F,DF,SDEL,3,1,2,4,MIN);
    S4=case8(S4,N,NM,R3,R13,R2,Z3,Z3L,Z13,Z13L,Z2,Z2L,F,DF,SDEL,3,1,4,2,MIN);
    S4=case8(S4,N,NM,R4,R14,R3,Z4,Z4L,Z14,Z14L,Z3,Z3L,F,DF,SDEL,4,1,2,3,MIN);
    S4=case8(S4,N,NM,R4,R14,R2,Z4,Z4L,Z14,Z14L,Z2,Z2L,F,DF,SDEL,4,1,3,2,MIN);
    S4=case8(S4,N,NM,R3,R23,R4,Z3,Z3L,Z23,Z23L,Z4,Z4L,F,DF,SDEL,3,2,1,4,MIN);
    S4=case8(S4,N,NM,R3,R23,R1,Z3,Z3L,Z23,Z23L,Z1,Z1L,F,DF,SDEL,3,2,4,1,MIN);
    S4=case8(S4,N,NM,R4,R24,R3,Z4,Z4L,Z24,Z24L,Z3,Z3L,F,DF,SDEL,4,2,1,3,MIN);
    S4=case8(S4,N,NM,R4,R24,R1,Z4,Z4L,Z24,Z24L,Z1,Z1L,F,DF,SDEL,4,2,3,1,MIN);
    S4=case8(S4,N,NM,R4,R34,R2,Z4,Z4L,Z34,Z34L,Z2,Z2L,F,DF,SDEL,4,3,1,2,MIN);
    S4=case8(S4,N,NM,R4,R34,R1,Z4,Z4L,Z34,Z34L,Z1,Z1L,F,DF,SDEL,4,3,2,1,MIN);
end
%}    
    
    S4=zeros(2,2,2,2);
    %S4split=zeros(8,2,2,2,2);
    orders = perms(1:4);
    for orderNum=1:24
        order=orders(orderNum,:);
        Q = [Q1,Q2,Q3,Q4];
        Qnew=[Q(:,order(1)),Q(:,order(2)),Q(:,order(3)),Q(:,order(4))];
        % Qnew is the reordered Q
        
        % Now calculate the eigenvalues 
        R1=-dot(Qnew(:,1),Qnew(:,1))/(2*d);
        R12=-dot(Qnew(:,1)+Qnew(:,2),Qnew(:,1)+Qnew(:,2))/(2*d);
        R4=-dot(Qnew(:,4),Qnew(:,4))/(2*d);
        
        
        Z1=exp(R1*NM);
        Z12=exp(R12*NM);
        Z4=exp(R4*NM);
        Z1L=Z1*LAM;
        Z12L=Z12*LAM;
        Z4L=Z4*LAM;
        
        S4=S4+case1(N,NM,R1,R12,R4,F);
        S4=S4+case2(N,NM,R1,R12,R4,Z1,Z1L,                F,order);
        S4=S4+case3(N,NM,R1,R12,R4,       Z12,Z12L,       F,order);
        S4=S4+case4(N,NM,R1,R12,R4,                Z4,Z4L,F,order);
        S4=S4+case5(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,       F,order);
        S4=S4+case6(N,NM,R1,R12,R4,       Z12,Z12L,Z4,Z4L,F,order);
        S4=S4+case7(N,NM,R1,R12,R4,Z1,Z1L,         Z4,Z4L,F,order);
        S4=S4+case8(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z4,Z4L,F,order);
        
%{        
        S4split(1,:,:,:,:)=squeeze(S4split(1,:,:,:,:))+case1(N,NM,R1,R12,R4,F);
        S4split(2,:,:,:,:)=squeeze(S4split(2,:,:,:,:))+case2(N,NM,R1,R12,R4,Z1,Z1L,                F,order);
        S4split(3,:,:,:,:)=squeeze(S4split(3,:,:,:,:))+case3(N,NM,R1,R12,R4,       Z12,Z12L,       F,order);
        S4split(4,:,:,:,:)=squeeze(S4split(4,:,:,:,:))+case4(N,NM,R1,R12,R4,                Z4,Z4L,F,order);
        S4split(5,:,:,:,:)=squeeze(S4split(5,:,:,:,:))+case5(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,       F,order);
        S4split(6,:,:,:,:)=squeeze(S4split(6,:,:,:,:))+case6(N,NM,R1,R12,R4,       Z12,Z12L,Z4,Z4L,F,order);
        S4split(7,:,:,:,:)=squeeze(S4split(7,:,:,:,:))+case7(N,NM,R1,R12,R4,Z1,Z1L,         Z4,Z4L,F,order);
        S4split(8,:,:,:,:)=squeeze(S4split(8,:,:,:,:))+case8(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z4,Z4L,F,order);
%}

    end
    
end




