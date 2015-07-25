function S4split=s4gaussian(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d)
%S4 is a four point correlation function
%For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA

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
    
    %new
    S4=zeros(2,2,2,2);
    S4split=zeros(8,2,2,2,2);
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
        
        
        S4split(1,:,:,:,:)=squeeze(S4split(1,:,:,:,:))+case1(N,NM,R1,R12,R4,F);
        S4split(2,:,:,:,:)=squeeze(S4split(2,:,:,:,:))+case2(N,NM,R1,R12,R4,Z1,Z1L,                F,order);
        S4split(3,:,:,:,:)=squeeze(S4split(3,:,:,:,:))+case3(N,NM,R1,R12,R4,       Z12,Z12L,       F,order);
        S4split(4,:,:,:,:)=squeeze(S4split(4,:,:,:,:))+case4(N,NM,R1,R12,R4,                Z4,Z4L,F,order);
        S4split(5,:,:,:,:)=squeeze(S4split(5,:,:,:,:))+case5(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,       F,order);
        S4split(6,:,:,:,:)=squeeze(S4split(6,:,:,:,:))+case6(N,NM,R1,R12,R4,       Z12,Z12L,Z4,Z4L,F,order);
        S4split(7,:,:,:,:)=squeeze(S4split(7,:,:,:,:))+case7(N,NM,R1,R12,R4,Z1,Z1L,         Z4,Z4L,F,order);
        S4split(8,:,:,:,:)=squeeze(S4split(8,:,:,:,:))+case8(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z4,Z4L,F,order);

    end
    
end

function S4=case1(N,NM,R1,R12,R3,F)
% Case 1: J1==J2==J3==J4
    E1=R1;
    E12=R12;
    E3=R3;
    
    % on same monomer
    offset=10^-12;
    if abs(E1)<offset
        E1=offset;
    end
    if abs(E12)<offset
        E12=offset;
    end
    if abs(E3)<offset
        E3=offset;
    end
    MIN=(10^-3)/NM;
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(NM*(E1+E12+E3)+5)/120;
    elseif max(abs([E1,E12]))<MIN
        valeq=chicken(E1,E12,E3,NM);
    elseif max(abs([E12,E3]))<MIN
        valeq=chicken(E12,E3,E1,NM);
    elseif max(abs([E1,E3]))<MIN
        valeq=chicken(E1,E3,E12,NM);
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
        valeq=(1/2).*E1.^(-4).*((-2).*(3+E1.*NM)+exp(1).^(E1.*NM).*(6+E1.*NM.*(( ...
         -4)+E1.*NM)));
    elseif (abs(E1-E12)<MIN )%&& abs(E12-E3)>=MIN)
        valeq=E1.^(-3).*(E1+(-1).*E3).^(-2).*E3.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
          NM)).*E3.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E3.^2.*NM+E1.^3.*((-1)+ ...
          exp(1).^(E3.*NM)+(-1).*E3.*NM)+(-1).*E1.*E3.^2.*((-3)+E3.*NM+exp( ...
          1).^(E1.*NM).*(3+E3.*NM)));
    elseif (abs(E1-E3)<MIN )%&& abs(E1-E12)>=MIN)
%         valeq=E1.^(-3).*(E1+(-1).*E12).^(-2).*E12.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
%           NM)).*E12.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E12.^2.*NM+E1.^3.*((-1)+ ...
%           exp(1).^(E12.*NM)+(-1).*E12.*NM)+(-1).*E1.*E12.^2.*((-3)+E12.*NM+ ...
%           exp(1).^(E1.*NM).*(3+E12.*NM)));
        valeq=((-3*E1*expl(4,NM*E1)+NM*E1^2*expl(3,NM*E1))*E12^2+...
          (2*expl(4,NM*E1)-NM*E1*expl(3,NM*E1))*E12^3+...
          E1^3*expl(4,NM*E12))/(E1^3*E12^2*(E1-E12)^2);
    elseif (abs(E12-E3)<MIN )%&& abs(E1-E3)>=MIN)
       valeq=E1.^(-2).*(E1+(-1).*E12).^(-2).*E12.^(-3).*(((-1)+exp(1).^(E1.*NM) ...
          ).*E12.^3+(-1).*E1.*E12.^3.*NM+E1.^2.*E12.*(3+2.*E12.*NM+exp(1).^( ...
          E12.*NM).*((-3)+E12.*NM))+(-1).*E1.^3.*(2+E12.*NM+exp(1).^(E12.* ...
          NM).*((-2)+E12.*NM)));
    else
       valeq=E1.^(-2).*(E1+(-1).*E12).^(-1).*E12.^(-2).*(E1+(-1).*E3).^(-1).*( ...
          E12+(-1).*E3).^(-1).*E3.^(-2).*(((-1)+exp(1).^(E1.*NM)).*E12.^2.*( ...
          E12+(-1).*E3).*E3.^2+E1.*E12.^2.*E3.^2.*((-1).*E12+E3).*NM+E1.^3.* ...
          ((-1).*((-1)+exp(1).^(E12.*NM)).*E3.^2+E12.*E3.^2.*NM+E12.^2.*(( ...
          -1)+exp(1).^(E3.*NM)+(-1).*E3.*NM))+E1.^2.*(((-1)+exp(1).^(E12.* ...
          NM)).*E3.^3+(-1).*E12.*E3.^3.*NM+E12.^3.*(1+(-1).*exp(1).^(E3.*NM) ...
          +E3.*NM)));
    end
    
    S4=zeros(2,2,2,2);
    S4(1,1,1,1)=S4(1,1,1,1)+F(1)*N*valeq;
    S4(2,2,2,2)=S4(2,2,2,2)+F(2)*N*valeq;
end

function S4=case2(N,NM,R1,R12,R3,Z1,Z1L,F,order)
% Case 2: J1 <J2==J3==J4
    if N<2
        S4=zeros(2,2,2,2);
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=Z1;
    ZE1L=Z1L;

    valeq=case4Int(E3,E12,E1,NM);
    % on different monomers
    valne1=twoSum(ZE1,N);
    valne2=twoSum(ZE1L,N);  
    
    valne=zeros(2,2,2);
    valne(1,:,:)=ones(2,2)*valne1;
    valne(2,:,:)=ones(2,2)*valne2;
    
    S4=BinomialSum(valeq,valne,order,F);
end

function S4=case3(N,NM,R1,R12,R3,Z12,Z12L,F,order)
% Case 3: J1==J2 <J3==J4
    if N<2
        S4=zeros(2,2,2,2);
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE12=Z12;
    ZE12L=Z12L;

    % on same monomer
    offset=10^-12;
    if abs(E1)<offset
        E1=offset;
    end
    if abs(E12)<offset
        E12=offset;
    end
    if abs(E3)<offset
        E3=offset;
    end
    MIN=(10^-4)/NM;
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(NM*(E1-E12+E3)+3)/12;
    elseif max(abs([E1,E12]))<MIN
        valeq=(-NM^2/(12*E3^3))*((NM^2*(2*E1-E12)+6*NM)*E3^2+...
                    2*NM*E3*((1+2*exp(NM*E3))*E12+(-expl(1,NM*E3))*E1)+...
                    6*(-expl(1,NM*E3))*(E3+E12));
    elseif max(abs([E12,E3]))<MIN
        valeq=(-NM^2/(12*E1^3))*((NM^2*(2*E3-E12)+6*NM)*E1^2+...
                    2*NM*E1*((1+2*exp(NM*E1))*E12+(-expl(1,NM*E1))*E3)+...
                    6*(-expl(1,NM*E1))*(E1+E12));
    elseif max(abs([E1,E3]))<MIN
        valeq=((2*NM^2+NM^3*(E1+E3))*exp(-NM*E12)*E12^3+...
                2*NM^2*(E1+E3)*E12^2+...
                (expl(1,-NM*E12))*(...
                      (3*NM^2*(E1+E3)+4*NM)*E12^2+...
                      4*NM*(E1+E3)*E12 ...
                      -2*(E1+E12+E3)*expl(1,NM*E12)))/(2*E12^5);
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(exp(1).^E1).^((-1).*NM).*log(exp(1).^E1).^(-4).*(1+(exp(1).^E1) ...
       .^NM.*((-1)+NM.*log(exp(1).^E1))).^2;
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*log(exp(1).^E1).^(-2).*(1+(exp(1).^E1) ...
       .^NM.*((-1)+NM.*log(exp(1).^E1))).*(((-1)+(exp(1).^E1).^NM).*log( ...
       exp(1).^E1).^(-1)+(1+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E3).^( ...
       -1)).*(log(exp(1).^E1)+(-1).*log(exp(1).^E3)).^(-1);
    elseif abs(E1-E3)<MIN
       valeq=(exp(1).^E12).^((-1).*NM).*log(exp(1).^E1).^(-2).*(log(exp(1).^E1) ...
       +(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^(-2).*(((-1)+( ...
       exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*((-1)+(exp(1).^E1).^NM).* ...
       log(exp(1).^E12)).^2;
    elseif abs(E12-E3)<MIN
       valeq=log(exp(1).^E1).^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^( ...
       -1).*log(exp(1).^E12).^(-3).*((-1).*((-1)+(exp(1).^E12).^NM).*log( ...
       exp(1).^E1)+((-1)+(exp(1).^E1).^NM).*log(exp(1).^E12)).*((-1)+( ...
       exp(1).^E12).^((-1).*NM)+NM.*log(exp(1).^E12));
    else
       valeq=(-1).*(exp(1).^E12).^((-1).*NM).*log(exp(1).^E1).^(-1).*(log(exp( ...
       1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1).^E12).^(-1).*((( ...
       -1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*((-1)+(exp(1).^E1) ...
       .^NM).*log(exp(1).^E12)).*(((-1)+(exp(1).^E12).^NM).*log(exp(1) ...
       .^E12).^(-1)+(1+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E3).^(-1)).*( ...
       log(exp(1).^E12)+(-1).*log(exp(1).^E3)).^(-1);
    end
    if valeq<0
        sprintf('E1=%g, E12=%g,E3=%g',E1,E12,E3)
        error('value cannot be negitive')
    end
    
    % on different monomers
    valne1=twoSum(ZE12,N);
    valne2=twoSum(ZE12L,N);
    
    valne=zeros(2,2,2);
    valne(:,1,:)=ones(2,2)*valne1;
    valne(:,2,:)=ones(2,2)*valne2;
    
    S4=BinomialSum(valeq,valne,order,F);
end

function S4=case4(N,NM,R1,R12,R3,Z3,Z3L,F,order)
    if N<2
        S4=zeros(2,2,2,2);
        return
    end
    % Case 4: J1==J2==J3 <J4
    E1=R1;
    E12=R12;
    E3=R3;
    ZE3=Z3;
    ZE3L=Z3L;
   
    valeq=case4Int(E1,E12,E3,NM);
    % on different monomers
    valne1=twoSum(ZE3,N);
    valne2=twoSum(ZE3L,N);

    valne=zeros(2,2,2);
    valne(:,:,1)=ones(2,2)*valne1;
    valne(:,:,2)=ones(2,2)*valne2;
    
    S4=BinomialSum(valeq,valne,order,F);
end

function S4=case5(N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,F,order)
    if N<3
        S4=zeros(2,2,2,2);
        return
    end
    % Case 5: J1 <J2 <J3==J4 
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];

    % on same monomer
    %{
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*NM.*log(exp(1) ...
       .^E1).^(-3).*(1+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1)));
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*NM.*log(exp(1) ...
       .^E1).^(-1).*(((-1)+(exp(1).^E1).^NM).*log(exp(1).^E1).^(-1)+(1+( ...
       -1).*(exp(1).^E3).^NM).*log(exp(1).^E3).^(-1)).*(log(exp(1).^E1)+( ...
       -1).*log(exp(1).^E3)).^(-1);
    elseif abs(E1-E3)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*(exp(1).^E12).^((-1).*NM).*((-1)+(exp(1) ...
       .^E1).^NM).*((exp(1).^E1).^NM+(-1).*(exp(1).^E12).^NM).*log(exp(1) ...
       .^E1).^(-1).*(((-1)+(exp(1).^E1).^NM).*log(exp(1).^E1).^(-1)+(1+( ...
       -1).*(exp(1).^E12).^NM).*log(exp(1).^E12).^(-1)).*(log(exp(1).^E1) ...
       +(-1).*log(exp(1).^E12)).^(-2);
    elseif abs(E12-E3)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*((exp(1).^E1) ...
       .^NM+(-1).*(exp(1).^E12).^NM).*log(exp(1).^E1).^(-1).*(log(exp(1) ...
       .^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1).^E12).^(-2).*((-1) ...
       +(exp(1).^E12).^((-1).*NM)+NM.*log(exp(1).^E12));
    else
       valeq=(exp(1).^E1).^((-1).*NM).*(exp(1).^E12).^((-1).*NM).*((-1)+(exp(1) ...
       .^E1).^NM).*((exp(1).^E1).^NM+(-1).*(exp(1).^E12).^NM).*log(exp(1) ...
       .^E1).^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*((( ...
       -1)+(exp(1).^E12).^NM).*log(exp(1).^E12).^(-1)+(1+(-1).*(exp(1) ...
       .^E3).^NM).*log(exp(1).^E3).^(-1)).*(log(exp(1).^E12)+(-1).*log( ...
       exp(1).^E3)).^(-1);
    end
    %}
    valeq=case6Int(E3,E12,E1,NM);
  
    % on different monomers
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT12=ZE12(I2);
            valne(I1,I2,1)=tripleSum(ZT1,ZT12,N);
            valne(I1,I2,2)=valne(I1,I2,1);
        end
    end
%     valeq=finite(valeq);
%     valne=finite(valne);
%{
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A1)*FV(A2)*FV(A3)*valne(1,1)+...
                          FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(1,2)+...
                          FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*valne(2,1)+...
                          FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(2,2));
                end
            end
        end
    end
%}
    
    S4=BinomialSum(valeq,valne,order,F);
    
end

function S4=case6(N,NM,R1,R12,R3,Z12,Z12L,Z3,Z3L,F,order)
    if N<3
        S4=zeros(2,2,2,2);
        return
    end
    % Case 6: J1==J2 <J3 <J4
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
   
    % on same monomer\
    valeq=case6Int(E1,E12,E3,NM);
    
    % on different monomers
    valne=zeros(2,2,2);
    for I2=1:2
        for I3=1:2
            ZT12=ZE12(I2);
            ZT3=ZE3(I3);
            valne(1,I2,I3)=tripleSum(ZT12,ZT3,N);
            valne(2,I2,I3)=valne(1,I2,I3);
        end
    end
%     valeq=finite(valeq);
%     valne=finite(valne);

    S4=BinomialSum(valeq,valne,order,F);
end

function S4=case7(N,NM,R1,R12,R3,Z1,Z1L,Z3,Z3L,F,order)
    if N<3
        S4=zeros(2,2,2,2);
        return
    end
    % Case 7: J1 <J2==J3 <J4
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE1=[Z1,Z1L];
    ZE3=[Z3,Z3L];

    % on same monomer
    MIN=(10^-4)/NM;
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(1/2).*(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).^2.* ...
       NM.^2.*log(exp(1).^E1).^(-2);
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*(1+(-1).*(exp( ...
       1).^E3).^((-1).*NM)).*log(exp(1).^E1).^(-1).*(log(exp(1).^E1)+(-1) ...
       .*log(exp(1).^E3)).^(-2).*log(exp(1).^E3).^(-1).*((exp(1).^E3) ...
       .^NM+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1)+(-1).*NM.*log( ...
       exp(1).^E3)));
    elseif abs(E1-E3)<MIN
       valeq=(exp(1).^E1).^((-2).*NM).*((-1)+(exp(1).^E1).^NM).^2.*log(exp(1) ...
       .^E1).^(-2).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*(( ...
       exp(1).^E12).^NM+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1)+(-1) ...
       .*NM.*log(exp(1).^E12)));
    elseif abs(E12-E3)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*(1+(-1).*(exp( ...
       1).^E12).^((-1).*NM)).*log(exp(1).^E1).^(-1).*(log(exp(1).^E1)+( ...
       -1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^(-1).*((exp(1) ...
       .^E1).^NM+(-1).*(exp(1).^E12).^NM+(-1).*(exp(1).^E12).^NM.*NM.* ...
       log(exp(1).^E1)+(exp(1).^E12).^NM.*NM.*log(exp(1).^E12));
    else
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*(1+(-1).*(exp( ...
       1).^E3).^((-1).*NM)).*log(exp(1).^E1).^(-1).*(log(exp(1).^E1)+(-1) ...
       .*log(exp(1).^E12)).^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E3)) ...
       .^(-1).*log(exp(1).^E3).^(-1).*((-1).*log(exp(1).^E12)+log(exp(1) ...
       .^E3)).^(-1).*(((exp(1).^E12).^NM+(-1).*(exp(1).^E3).^NM).*log( ...
       exp(1).^E1)+((-1).*(exp(1).^E1).^NM+(exp(1).^E3).^NM).*log(exp(1) ...
       .^E12)+((exp(1).^E1).^NM+(-1).*(exp(1).^E12).^NM).*log(exp(1).^E3) ...
       );
    end

%     valeq=expl(1,-NM*E1)*expl(1,-NM*E3)*(...
%           (expl(2,NM*E3)-expl(2,NM*E12))*(E1-E3)+...
%           (expl(2,NM*E1)-expl(2,NM*E3))*(E12-E3))/...
%           ( E1*E3*(E12-E3)*(E1-E3)*(E1-E12) );
    
    % on different monomers
    
    valne=zeros(2,2,2);
    for I1=1:2
        for I3=1:2
            ZT1=ZE1(I1);
            ZT3=ZE3(I3);
            valne(I1,1,I3)=tripleSum(ZT1,ZT3,N);
            valne(I1,2,I3)=valne(I1,1,I3);
        end
    end

    S4=BinomialSum(valeq,valne,order,F);

end

function S4=case8(N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,order)
    if N<4
        S4=zeros(2,2,2,2);
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
    MIN=(10^-4)/NM;
    if max(abs([E1,E12,E3]))<MIN
       valeq=NM^4;
    elseif max(abs([E1,E12]))<MIN
        valeq=(NM^2/(2*E3^3))*(2*(E12+E3)*(expl(2,NM*E3)+expl(2,-NM*E3))...
            -2*NM*E12*E3*sinh(NM*E3));
    elseif max(abs([E12,E3]))<MIN
        valeq=(NM^2/(2*E1^3))*(2*(E12+E1)*(expl(2,NM*E1)+expl(2,-NM*E1))...
            -2*NM*E12*E1*sinh(NM*E1));
    elseif max(abs([E1,E3]))<MIN
        valeq=(NM^2/(E12^3))*(4*(E1+E12+E3)*sinh(0.5*NM*E12)^2 ...
            -NM*(E1+E3)*E12*sinh(NM*E12));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=2.*E1.^(-2).*NM.^2.*((-1)+cosh(E1.*NM));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1) ...
          .*(E1+(-1).*E3).^(-1).*E3.^(-1).*NM;
    elseif (abs(E1-E3)<MIN && abs(E1-E12)>=MIN)
       valeq=16.*E1.^(-2).*(E1+(-1).*E12).^(-2).*sinh((1/2).*E1.*NM).^2.*sinh(( ...
        1/2).*(E1+(-1).*E12).*NM).^2;
    elseif (abs(E12-E3)<MIN && abs(E1-E3)>=MIN)
       valeq=exp(1).^((-1).*(E1+E12).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1).^( ...
          E1.*NM)+(-1).*exp(1).^(E12.*NM)).*((-1)+exp(1).^(E12.*NM)).*E1.^( ...
          -1).*(E1+(-1).*E12).^(-1).*E12.^(-1).*NM;
    else
%        valeq=exp(1).^((-1).*(E1+E12+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1) ...
%           .^(E1.*NM)+(-1).*exp(1).^(E12.*NM)).*(exp(1).^(E12.*NM)+(-1).*exp( ...
%           1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1).*(E1+(-1).*E12) ...
%           .^(-1).*(E12+(-1).*E3).^(-1).*E3.^(-1);    
        valeq=2*( -coshl(4,NM*E1)+coshl(4,NM*E12)-coshl(4,NM*E3)-coshl(4,NM*(E1-E12))+...
           coshl(4,NM*(E1-E3))-coshl(4,NM*(E12-E3))+coshl(4,NM*(-E1+E12-E3)) )...
           *( (E12-E3)*E3*E1*(E12-E1) )^(-1);    
    end
    
    % on different monomers
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            for I3=1:2
                ZT1=ZE1(I1);
                ZT12=ZE12(I2);
                ZT3=ZE3(I3);
                valne(I1,I2,I3)=case8JPart(ZT1,ZT12,ZT3,N);
            end
        end
    end
    
    S4=BinomialSum(valeq,valne,order,F);
    
end

% vvvvvvv helper functions vvvvvvvv

function out=chicken(a,b,c,NM)
  out=(-NM^3*(NM*(a+b)+4)*c^4 ...
     -4*NM^2*(NM*(a+b)+3)*c^3 ...
    -12*NM^1*(NM*(a+b)+2)*c^2 ...
    +24*(expl(1,NM*c)-NM*(a+b))*c ...
    +24*expl(1,NM*c)*(a+b))/(24*c^5);
end

function valeq=case6Int(E1,E12,E3,NM)
    offset=10^-12;
    if abs(E1)<offset
        E1=offset;
    end
    if abs(E12)<offset
        E12=offset;
    end
    if abs(E3)<offset
        E3=offset;
    end
    MIN=(10^-4)/NM;
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(2*NM*E1-NM*E12+6)/12;
        flag=1;
    elseif max(abs([E1,E12]))<MIN
        valeq=NM^2*exp(NM*E3)*expl(1,-NM*E3)*(3*(E12+E3)*expl(2,-NM*E3)-3*NM*E3^2 ...
              +NM*(E1+E12)*E3*expl(1,-NM*E3))/(6*E3^3);
        flag=2;
    elseif max(abs([E12,E3]))<MIN
        valeq=(NM^2/(E1^3))*((E1+E12)*expl(2,NM*E1)-0.5*NM*E1*E12*expl(1,NM*E1));
        flag=3;
    elseif max(abs([E1,E3]))<MIN
        valeq=(NM/(2*E12^4))*(2*(E1+E12+E3)*...
               (expl(3,NM*E12)+expl(3,-NM*E12)+NM*E12*expl(2,-NM*E12))+...
               NM^2*(E1+E3)*E12^2*expl(1,-NM*E12)+...
               NM*E12*E3*(expl(2,-NM*E12)-expl(2,NM*E12)));
        flag=4;
    elseif(abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*NM.*log(exp(1) ...
       .^E1).^(-3).*(1+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1)));
       flag=5;
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*(1+(-1).*(exp(1).^E3).^((-1).*NM)).*(( ...
       exp(1).^E1).^NM+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E1).^(-2).*( ...
       1+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1))).*(log(exp(1).^E1)+ ...
       (-1).*log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1);
       flag=6;
    elseif abs(E1-E3)<MIN
       valeq=(exp(1).^E12).^((-1).*NM).*(1+(-1).*(exp(1).^E1).^((-1).*NM)).*(( ...
       exp(1).^E1).^NM+(-1).*(exp(1).^E12).^NM).*log(exp(1).^E1).^(-2).*( ...
       log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^( ...
       -1).*((-1).*((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+((-1)+(exp( ...
       1).^E1).^NM).*log(exp(1).^E12));
       flag=7;
    elseif abs(E12-E3)<MIN
       valeq=(1+(-1).*(exp(1).^E12).^((-1).*NM)).*NM.*log(exp(1).^E1).^(-1).*( ...
       log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1).^E12).^( ...
       -2).*((-1).*((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+((-1)+(exp( ...
       1).^E1).^NM).*log(exp(1).^E12));
       flag=8;
    else
        valeq=expl(1,NM*E1)*(expl(1,-NM*E1)-expl(1,-NM*E12))*...
              ( expl(2,NM*E3)*E12 - expl(2,NM*E12)*E3 )/...
              (E1*E12*E3*(E12-E3)*(E1-E12));
        flag=9;
%        valeq=(-1).*(exp(1).^E12).^((-1).*NM).*(1+(-1).*(exp(1).^E3).^((-1).*NM) ...
%        ).*((exp(1).^E12).^NM+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E1).^( ...
%        -1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1) ...
%        .^E12).^(-1).*(((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*(( ...
%        -1)+(exp(1).^E1).^NM).*log(exp(1).^E12)).*(log(exp(1).^E12)+(-1).* ...
%        log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1);
    end
    if valeq<0
        sprintf('a=%g, b=%g, c=%g, valeq=%g, flag=%d',E1,E12,E3,valeq,flag)
        error('valeq<0')
    end
end

function valeq=case4Int(E1,E12,E3,NM)
% This function calculates the intigral for case 4
% This is the same intigral as case 2 with arguements exchanged
    offset=10^-12;
    if abs(E1)<offset
        E1=offset;
    end
    if abs(E12)<offset
        E12=offset;
    end
    if abs(E3)<offset
        E3=offset;
    end
    % on same monomer
    MIN=(10^-4)/NM;
    if max(abs([E1,E12,E3]))<MIN  % 0~E1~E12~E3
        valeq=NM^4*(NM*(E1+E12-E3)+4)/24;
    elseif max(abs([E1,E3]))<MIN   % 0~E1~E3
        valeq=((NM^4*(E3-2*E1)-6*NM^3)*E12^3+...
            (-6*E1*NM^3-12*NM^2)*E12^2+...
            ((12*NM-6*NM^2*E3)*exp(NM*E12)-6*NM^2 *(E3+2*E1)-12*NM)*E12+...
            12*NM*(E1+E3)*(exp(NM*E12)-1))/(12*E12^4);
    elseif max(abs([E1,E12]))<MIN  % 0~E1~E12
        valeq=(expl(1,-NM*E3)*((E1+E12)*((NM*E3)^3+3*(NM*E3)^2+6*NM*E3)...
            +3*NM*E3^2*(2+NM*E3))...
            +6*(expl(2,NM*E3)+expl(2,-NM*E3))*(E1+E12+E3))/(6*E3^5);
    elseif max(abs([E12,E3]))<MIN % 0~E12~E3
        valeq=-NM*(  (NM*(2*E12-E3)+6)*NM^2*E1^3+...
           6*(E12*NM^2+2*NM)*E1^2+...
           6*(NM*(2*E12+E3)+2+(NM*E3-2)*exp(NM*E1))*E1+12*(E3+E12)*(1-exp(NM*E1)) )/...
           (12*E1^4);
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(1/2).*(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*log(exp( ...
       1).^E1).^(-4).*((-2)+(exp(1).^E1).^NM.*(2+NM.*log(exp(1).^E1).*(( ...
       -2)+NM.*log(exp(1).^E1))));
    elseif abs(E1-E12)<MIN
       valeq=(1+(-1).*(exp(1).^E3).^((-1).*NM)).*log(exp(1).^E1).^(-2).*(log( ...
       exp(1).^E1)+(-1).*log(exp(1).^E3)).^(-2).*log(exp(1).^E3).^(-2).*( ...
       ((-1)+(exp(1).^E1).^NM).*log(exp(1).^E3).^2+(-1).*log(exp(1).^E1) ...
       .*log(exp(1).^E3).*((-2)+2.*(exp(1).^E1).^NM+(exp(1).^E1).^NM.* ...
       NM.*log(exp(1).^E3))+log(exp(1).^E1).^2.*((-1)+(exp(1).^E3).^NM+( ...
       exp(1).^E1).^NM.*NM.*log(exp(1).^E3)));
    elseif abs(E1-E3)<MIN
       valeq=(1+(-1).*(exp(1).^E1).^((-1).*NM)).*log(exp(1).^E1).^(-3).*(log( ...
       exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^(-1) ...
       .*(((-1)+(exp(1).^E1).^NM).*log(exp(1).^E12).^2+(-1).*log(exp(1) ...
       .^E1).*log(exp(1).^E12).*((-2)+2.*(exp(1).^E1).^NM+(exp(1).^E1) ...
       .^NM.*NM.*log(exp(1).^E12))+log(exp(1).^E1).^2.*((-1)+(exp(1) ...
       .^E12).^NM+(exp(1).^E1).^NM.*NM.*log(exp(1).^E12)));
    elseif abs(E12-E3)<MIN
       valeq=(1+(-1).*(exp(1).^E12).^((-1).*NM)).*log(exp(1).^E1).^(-1).*(log( ...
       exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^(-3) ...
       .*(((-1)+(exp(1).^E1).^NM).*log(exp(1).^E12).^2+log(exp(1).^E1) ...
       .^2.*((-1)+(exp(1).^E12).^NM+(-1).*(exp(1).^E12).^NM.*NM.*log(exp( ...
       1).^E12))+log(exp(1).^E1).*log(exp(1).^E12).*(2+(-2).*(exp(1) ...
       .^E12).^NM+(exp(1).^E12).^NM.*NM.*log(exp(1).^E12)));
    else
       valeq=(1+(-1).*(exp(1).^E3).^((-1).*NM)).*log(exp(1).^E1).^(-1).*(log( ...
       exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1).^E12).^(-1) ...
       .*(log(exp(1).^E1)+(-1).*log(exp(1).^E3)).^(-1).*(log(exp(1).^E12) ...
       +(-1).*log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-2).*(((-1)+(exp( ...
       1).^E1).^NM).*log(exp(1).^E12).*(log(exp(1).^E12)+(-1).*log(exp(1) ...
       .^E3)).*log(exp(1).^E3)+log(exp(1).^E1).^2.*(((-1)+(exp(1).^E3) ...
       .^NM).*log(exp(1).^E12)+(-1).*((-1)+(exp(1).^E12).^NM).*log(exp(1) ...
       .^E3))+log(exp(1).^E1).*((-1).*((-1)+(exp(1).^E3).^NM).*log(exp(1) ...
       .^E12).^2+((-1)+(exp(1).^E12).^NM).*log(exp(1).^E3).^2));     
%         x= expl(3,NM*E1 )*E3 -expl(3,NM*E3 )*E1 ;
%         y= expl(3,NM*E3 )*E12-expl(3,NM*E12)*E3 ;
%         z= expl(3,NM*E12)*E1 -expl(3,NM*E1 )*E12; 
%         valeq=expl(1,-NM*E3)*(x*E12^2+y*E1^2+z*E3^2)/...
%           (-E1*E12*E3^2*(E12-E3)*(E1-E3)*(E1-E12));
    end 
    if valeq<0
        sprintf('E1=%g, E12=%g, E3=%g, valeq=%g',E1,E12,E3,valeq)
        error('value should not be negitive')
    end
end

function out=twoSum(a,N)
   if N<2
       out=0;
       return
   elseif a<0.8
       out=a*(N-N*a+a^N-1)/((1-a)^2);
   elseif abs(1-a)<10^-13
       out=N*(N+1)/2 -N;
   else
       out=a*(-N*expl(2,log(a))+expl(2,N*log(a)))/((a-1)^2);
   end
end

function xf=finite(x)
    MIN=-1e5;
    MAX=1e5;
    xf=x;
    xf(xf>MAX | xf<MIN)=0;
    xf(isnan(xf))=0;
    error('you shouldn`t use this')
end

function dif=BinomialSum(valeq,valne,order,F)
    dif=zeros(2,2,2,2);
    for a1=1:2
        for a2=1:2
            for a3=1:2
                for a4=1:2
                    aspace=[a1,a2,a3,a4];
                    apath=aspace(order); %reorder alphas into new order
                    FV=[F(apath(1)),F(apath(2)),F(apath(3)),F(apath(4))];
                    Delta=[apath(1)==apath(2),apath(2)==apath(3),apath(3)==apath(4)]*2-1;  % +1 is same, -1 if different
                    
                    bin1=[FV(2),Delta(1)*(1-FV(1))]; % bin stanods for binomial
                    bin2=[FV(3),Delta(2)*(1-FV(2))];
                    bin3=[FV(4),Delta(3)*(1-FV(3))];
                    
                    dif=dif+valeq*FV(1)*(...
                        bin1(1)*bin2(1)*bin3(1)*valne(1,1,1)+...
                        bin1(1)*bin2(1)*bin3(2)*valne(1,1,2)+...
                        bin1(1)*bin2(2)*bin3(1)*valne(1,2,1)+...
                        bin1(1)*bin2(2)*bin3(2)*valne(1,2,2)+...
                        bin1(2)*bin2(1)*bin3(1)*valne(2,1,1)+...
                        bin1(2)*bin2(1)*bin3(2)*valne(2,1,2)+... 
                        bin1(2)*bin2(2)*bin3(1)*valne(2,2,1)+...
                        bin1(2)*bin2(2)*bin3(2)*valne(2,2,2));
%{                    
%                     dif2=SDEL(a1,a2,a3,a4)*...
%                     valeq*(FV(A1)*FV(A2)             *FV(A3)             *FV(A4)             *valne(1,1,1)+...
%                            FV(A1)*FV(A2)             *FV(A3)             *Delta(3)*(1-FV(A3))*valne(1,1,2)+...
%                            FV(A1)*FV(A2)             *Delta(2)*(1-FV(A2))*FV(A4)             *valne(1,2,1)+...
%                            FV(A1)*FV(A2)             *Delta(2)*(1-FV(A2))*Delta(3)*(1-FV(A3))*valne(1,2,2)+...
%                            FV(A1)*Delta(1)*(1-FV(A1))*FV(A3)             *FV(A4)             *valne(2,1,1)+...
%                            FV(A1)*Delta(1)*(1-FV(A1))*FV(A3)             *Delta(3)*(1-FV(A3))*valne(2,1,2)+...
%                            FV(A1)*Delta(1)*(1-FV(A1))*Delta(2)*(1-FV(A2))*FV(A4)             *valne(2,2,1)+...
%                            FV(A1)*Delta(1)*(1-FV(A1))*Delta(2)*(1-FV(A2))*Delta(3)*(1-FV(A3))*valne(2,2,2));

                    DF=[1,-1];
                    FV=[F(a1),F(a2),F(a3),F(a4)];
                    DFV=[DF(a1),DF(a2),DF(a3),DF(a4)];
                    
                    A1=order(1);A2=order(2);A3=order(3);A4=order(4);
                    
                    difOld=1*...  used to have SDEL(a1,a2,a3,a4)
                      valeq*(FV(A1)*FV(A2)*FV(A3)*FV(A4)*valne(1,1,1)+...
                           FV(A1)*FV(A2)*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,1,2)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(1,2,1)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,2,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*FV(A4)*valne(2,1,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,1,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(2,2,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,2,2));
                    
                    if abs((dif - difOld)/difOld)>10^-13  && ~isnan(difOld)
                        
                        format long
                        FV(A1)* bin1(1)*bin2(1)*bin3(2)
                        FV(A1)*FV(A2)*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))
                        
                        dif
                        difOld
                        order
                        [a1,a2,a3,a4]
                        error('dif ~= difOld')
                    end
%}                        
                    %S4(a1,a2,a3,a4)=S4(a1,a2,a3,a4)+dif;
                end
            end
        end
    end
end


