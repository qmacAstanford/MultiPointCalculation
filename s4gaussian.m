function S4=s4gaussian(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d)
%S4 is a four point correlation function
%For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA

S4=zeros(2,2,2,2);
MIN=1e-6;
% MIN=1e-10;
digits(5);

% Reset Qs to column vectors if entered as rows
if isrow(Q1)==1
    Q1=transpose(Q1);
    Q2=transpose(Q2);
    Q3=transpose(Q3);
    Q4=transpose(Q4);
end

% Begin calculation of s4
if sum(power(Q1+Q2+Q3+Q4,2)) <= MIN

    % Evaluate the quantities for s4 calculation
    FB=1-FA;
    F=[FA,FB];
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
    CasesIncluded=[1,1,1,1,1,1,1,1];
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
        S4=case6(S4,N,NM,R2,R24,R3,Z24,Z24L,Z3,Z3L,F,DF,SDEL,4,1,3,MIN);
        S4=case6(S4,N,NM,R4,R24,R3,Z24,Z24L,Z3,Z3L,F,DF,SDEL,2,1,3,MIN);
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
        S4=case7(S4,N,NM,R2,R12,R4,Z2,Z2L,Z4,Z4L,F,DF,SDEL,2,1,4,MIN);
        S4=case7(S4,N,NM,R2,R13,R4,Z2,Z2L,Z4,Z4L,F,DF,SDEL,2,3,4,MIN);
        % J4 <J1==J3 <J2
        S4=case7(S4,N,NM,R4,R12,R2,Z4,Z4L,Z2,Z2L,F,DF,SDEL,4,1,2,MIN);
        S4=case7(S4,N,NM,R4,R13,R2,Z4,Z4L,Z2,Z2L,F,DF,SDEL,4,3,2,MIN);
        
        % J2 <J1==J4 <J3 (two terms)
        SDEL=zeros(2,2,2,2);
        SDEL(1,1,1,1)=1;
        SDEL(2,2,2,2)=1;
        SDEL(1,1,2,1)=1;
        SDEL(1,2,1,1)=1;
        SDEL(1,2,2,1)=1;
        SDEL(2,1,2,2)=1;
        SDEL(2,2,1,2)=1;
        SDEL(2,1,1,2)=1;
        S4=case7(S4,N,NM,R2,R12,R3,Z2,Z2L,Z3,Z3L,F,DF,SDEL,2,1,3,MIN);
        S4=case7(S4,N,NM,R2,R24,R3,Z2,Z2L,Z3,Z3L,F,DF,SDEL,2,4,3,MIN);
        % J3 <J1==J4 <J2
        S4=case7(S4,N,NM,R3,R13,R2,Z3,Z3L,Z2,Z2L,F,DF,SDEL,3,1,2,MIN);
        S4=case7(S4,N,NM,R3,R34,R2,Z3,Z3L,Z2,Z2L,F,DF,SDEL,3,4,2,MIN);
        
        % J1 <J3==J4 <J2 (two terms)
        SDEL=zeros(2,2,2,2);
        SDEL(1,1,1,1)=1;
        SDEL(2,2,2,2)=1;
        SDEL(2,1,1,1)=1;
        SDEL(1,2,1,1)=1;
        SDEL(2,2,1,1)=1;
        SDEL(1,2,2,2)=1;
        SDEL(1,2,2,2)=1;
        SDEL(1,1,2,2)=1;
        S4=case7(S4,N,NM,R1,R13,R2,Z1,Z1L,Z2,Z2L,F,DF,SDEL,1,3,2,MIN);
        S4=case7(S4,N,NM,R1,R14,R2,Z1,Z1L,Z2,Z2L,F,DF,SDEL,1,4,2,MIN);
        % J2 <J3==J4 <J1
        S4=case7(S4,N,NM,R2,R23,R1,Z2,Z2L,Z1,Z1L,F,DF,SDEL,2,3,1,MIN);
        S4=case7(S4,N,NM,R2,R24,R1,Z2,Z2L,Z1,Z1L,F,DF,SDEL,2,4,1,MIN);
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
end
end

function S4=case1(S4,N,NM,R1,R12,R3,F,MIN)
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
    MIN=(10^-2)/NM;
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
    %valeq=finite(valeq);
    
    S4(1,1,1,1)=S4(1,1,1,1)+2*F(1)*N*valeq;
    S4(2,2,2,2)=S4(2,2,2,2)+2*F(2)*N*valeq;
end
function out=chicken(a,b,c,NM)
  out=(-NM^3*(NM*(a+b)+4)*c^4 ...
     -4*NM^2*(NM*(a+b)+3)*c^3 ...
    -12*NM^1*(NM*(a+b)+2)*c^2 ...
    +24*(expl(1,NM*c)-NM*(a+b))*c ...
    +24*expl(1,NM*c)*(a+b))/(24*c^5);
end
function S4=case2(S4,N,NM,R1,R12,R3,Z1,Z1L,F,DF,SDEL,A1,A2,MIN)
% Case 2: J1 <J2==J3==J4
    if N<2
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=Z1;
    ZE1L=Z1L;
    %{
    % on same monomer
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(NM*(E1-E12-E3)+4)/24;
    elseif (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(1/2).*(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*log(exp( ...
       1).^E1).^(-4).*((-2)+(exp(1).^E1).^NM.*(2+NM.*log(exp(1).^E1).*(( ...
       -2)+NM.*log(exp(1).^E1))));
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*log(exp(1).^E1) ...
       .^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E3)).^(-2).*(((-1)+( ...
       exp(1).^E1).^NM).*log(exp(1).^E1).^(-2).*log(exp(1).^E3)+log(exp( ...
       1).^E1).^(-1).*(2+(-2).*(exp(1).^E1).^NM+(-1).*(exp(1).^E1).^NM.* ...
       NM.*log(exp(1).^E3))+log(exp(1).^E3).^(-1).*((-1)+(exp(1).^E3) ...
       .^NM+(exp(1).^E1).^NM.*NM.*log(exp(1).^E3)));
    elseif abs(E1-E3)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*log(exp(1).^E1) ...
       .^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*(((-1)+( ...
       exp(1).^E1).^NM).*log(exp(1).^E1).^(-2).*log(exp(1).^E12)+log(exp( ...
       1).^E1).^(-1).*(2+(-2).*(exp(1).^E1).^NM+(-1).*(exp(1).^E1).^NM.* ...
       NM.*log(exp(1).^E12))+log(exp(1).^E12).^(-1).*((-1)+(exp(1).^E12) ...
       .^NM+(exp(1).^E1).^NM.*NM.*log(exp(1).^E12)));
    elseif abs(E12-E3)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*log(exp(1).^E1) ...
       .^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*(((-1)+( ...
       exp(1).^E1).^NM).*log(exp(1).^E1).^(-1)+log(exp(1).^E12).^(-2).*( ...
       log(exp(1).^E1).*((-1)+(exp(1).^E12).^NM+(-1).*(exp(1).^E12).^NM.* ...
       NM.*log(exp(1).^E12))+log(exp(1).^E12).*(2+(-2).*(exp(1).^E12) ...
       .^NM+(exp(1).^E12).^NM.*NM.*log(exp(1).^E12))));
    else
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*log(exp(1).^E1) ...
       .^(-2).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1) ...
       .^E12).^(-1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E3)).^(-1).*(log( ...
       exp(1).^E12)+(-1).*log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1).* ...
       (((-1)+(exp(1).^E1).^NM).*log(exp(1).^E12).*(log(exp(1).^E12)+(-1) ...
       .*log(exp(1).^E3)).*log(exp(1).^E3)+log(exp(1).^E1).^2.*(((-1)+( ...
       exp(1).^E3).^NM).*log(exp(1).^E12)+(-1).*((-1)+(exp(1).^E12).^NM) ...
       .*log(exp(1).^E3))+log(exp(1).^E1).*((-1).*((-1)+(exp(1).^E3).^NM) ...
       .*log(exp(1).^E12).^2+((-1)+(exp(1).^E12).^NM).*log(exp(1).^E3) ...
       .^2));
    end
    %}
    valeq=case4Int(E3,E12,E1,NM);
    % on different monomers
%     valne1=(-1).*((-1)+ZE1).^(-2).*ZE1.*(1+N.*((-1)+ZE1)+(-1).*ZE1.^N);
%     valne2=(-1).*((-1)+ZE1L).^(-2).*ZE1L.*(1+N.*((-1)+ZE1L)+(-1).*ZE1L.^N);
    valne1=twoSum(ZE1,N);
    valne2=twoSum(ZE1L,N);  
    
%     valeq=finite(valeq);
%     valne1=finite(valne1);
%     valne2=finite(valne2);
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A1)*FV(A2)*valne1+...
                         FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*valne2);
                end
            end
        end
    end
end

function S4=case3(S4,N,NM,R1,R12,R3,Z12,Z12L,F,DF,SDEL,A2,A3,MIN)
% Case 3: J1==J2 <J3==J4
    if N<2
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
    MIN=10^-5;
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
    %valne1=(-1).*((-1)+ZE12).^(-2).*ZE12.*(1+N.*((-1)+ZE12)+(-1).*ZE12.^N);
    valne1=twoSum(ZE12,N);
    %valne2=(-1).*((-1)+ZE12L).^(-2).*ZE12L.*(1+N.*((-1)+ZE12L)+(-1).*ZE12L.^N);
    valne2=twoSum(ZE12L,N);
    
%     valeq=finite(valeq);
%    valne1=finite(valne1);
%    valne2=finite(valne2);
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A2)*FV(A3)*valne1+...
                         FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne2);
                end
            end
        end
    end
end

function S4=case4(S4,N,NM,R1,R12,R3,Z3,Z3L,F,DF,SDEL,A3,A4,MIN)
    if N<2
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
    %valne1=(-1).*((-1)+ZE3).^(-2).*ZE3.*(1+N.*((-1)+ZE3)+(-1).*ZE3.^N)
    valne1=twoSum(ZE3,N);
    %valne2=(-1).*((-1)+ZE3L).^(-2).*ZE3L.*(1+N.*((-1)+ZE3L)+(-1).*ZE3L.^N);
    valne2=twoSum(ZE3L,N);

    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A3)*FV(A4)*valne1+...
                         FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne2);
                end
            end
        end
    end
end

function valeq=case4Int(E1,E12,E3,NM)
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

function S4=case5(S4,N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,F,DF,SDEL,A1,A2,A3,MIN)
    if N<3
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
    MIN=5e-3;
    valne=zeros(2,2);
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT12=ZE12(I2);
            valne(I1,I2)=tripleSum(ZT1,ZT12,N);
%{            
            if abs(ZT1-ZT12)<MIN
                valne(I1,I2)=((-1)+ZT1).^(-3).*ZT1.*(2.*ZT1+(-1).*N.*ZT1+N.*ZT1.^2+(-1).*N.* ...
                  ZT1.^N+(-2).*ZT1.^(1+N)+N.*ZT1.^(1+N));
            else
                valne(I1,I2)=((-1)+ZT1).^(-2).*ZT1.*(ZT1+(-1).*ZT12).^(-1).*((-1)+ZT12).^(-2).* ...
                  ZT12.*((-2).*ZT1+N.*ZT1+ZT1.^2+(-1).*N.*ZT1.^2+ZT1.^N+2.*ZT12+(-1) ...
                  .*N.*ZT12+N.*ZT1.^2.*ZT12+(-2).*ZT1.^N.*ZT12+(-1).*ZT12.^2+N.* ...
                  ZT12.^2+(-1).*N.*ZT1.*ZT12.^2+ZT1.^N.*ZT12.^2+(-1).*ZT12.^N+2.* ...
                  ZT1.*ZT12.^N+(-1).*ZT1.^2.*ZT12.^N);
            end
%}            
        end
    end
%     valeq=finite(valeq);
%     valne=finite(valne);

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
end

function S4=case6(S4,N,NM,R1,R12,R3,Z12,Z12L,Z3,Z3L,F,DF,SDEL,A2,A3,A4,MIN)
    if N<3
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
    MIN=5e-3;
    valne=zeros(2,2);
    for I1=1:2
        for I2=1:2
            ZT12=ZE12(I1);
            ZT3=ZE3(I2);
            valne(I1,I2)=tripleSum(ZT12,ZT3,N);
        end
    end
%     valeq=finite(valeq);
%     valne=finite(valne);

    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A2)*FV(A3)*FV(A4)*valne(1,1)+...
                          FV(A2)*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,2)+...
                          FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(2,1)+...
                          FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,2));
                end
            end
        end
    end
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
    
    MIN=10^-5;
    if max(abs([E1,E12,E3]))<MIN
        valeq=NM^4*(2*NM*E1-NM*E12+6)/12;
    elseif max(abs([E1,E12]))<MIN
        valeq=NM^2*exp(NM*E3)*expl(1,-NM*E3)*(3*(E12+E3)*expl(2,-NM*E3)-3*NM*E3^2 ...
              +NM*(E1+E12)*E3*expl(1,-NM*E3))/(6*E3^3);
    elseif max(abs([E12,E3]))<MIN
        valeq=(NM^2/(E1^3))*((E1+E12)*expl(2,NM*E1)-0.5*NM*E1*E12*expl(1,NM*E1));
    elseif max(abs([E1,E3]))
        valeq=(NM/(2*E12^4))*(2*(E1+E12+E3)*...
               (expl(3,NM*E12)+expl(3,-NM*E12)+NM*E12*expl(2,-NM*E12))+...
               NM^2*(E1+E3)*E12^2*expl(1,-NM*E12)+...
               NM*E12*E3*(expl(2,-NM*E12)-expl(2,NM*E12)));
    elseif(abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(exp(1).^E1).^((-1).*NM).*((-1)+(exp(1).^E1).^NM).*NM.*log(exp(1) ...
       .^E1).^(-3).*(1+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1)));
    elseif abs(E1-E12)<MIN
       valeq=(exp(1).^E1).^((-1).*NM).*(1+(-1).*(exp(1).^E3).^((-1).*NM)).*(( ...
       exp(1).^E1).^NM+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E1).^(-2).*( ...
       1+(exp(1).^E1).^NM.*((-1)+NM.*log(exp(1).^E1))).*(log(exp(1).^E1)+ ...
       (-1).*log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1);
    elseif abs(E1-E3)<MIN
       valeq=(exp(1).^E12).^((-1).*NM).*(1+(-1).*(exp(1).^E1).^((-1).*NM)).*(( ...
       exp(1).^E1).^NM+(-1).*(exp(1).^E12).^NM).*log(exp(1).^E1).^(-2).*( ...
       log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-2).*log(exp(1).^E12).^( ...
       -1).*((-1).*((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+((-1)+(exp( ...
       1).^E1).^NM).*log(exp(1).^E12));
    elseif abs(E12-E3)<MIN
       valeq=(1+(-1).*(exp(1).^E12).^((-1).*NM)).*NM.*log(exp(1).^E1).^(-1).*( ...
       log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1).^E12).^( ...
       -2).*((-1).*((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+((-1)+(exp( ...
       1).^E1).^NM).*log(exp(1).^E12));
    else
        valeq=expl(1,NM*E1)*(expl(1,-NM*E1)-expl(1,-NM*E12))*...
              ( expl(2,NM*E3)*E12 - expl(2,NM*E12)*E3 )/...
              (E1*E12*E3*(E12-E3)*(E1-E12));
%        valeq=(-1).*(exp(1).^E12).^((-1).*NM).*(1+(-1).*(exp(1).^E3).^((-1).*NM) ...
%        ).*((exp(1).^E12).^NM+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E1).^( ...
%        -1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1) ...
%        .^E12).^(-1).*(((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*(( ...
%        -1)+(exp(1).^E1).^NM).*log(exp(1).^E12)).*(log(exp(1).^E12)+(-1).* ...
%        log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1);
    end
end

function S4=case7(S4,N,NM,R1,R12,R3,Z1,Z1L,Z3,Z3L,F,DF,SDEL,A1,A2,A4,MIN)
    if N<3
        return
    end
    % Case 7: J1 <J2==J3 <J4
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE1=[Z1,Z1L];
    ZE3=[Z3,Z3L];

    % on same monomer
    MIN=10^-5;
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
    
    valne=zeros(2,2);
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT3=ZE3(I2);
            valne(I1,I2)=tripleSum(ZT1,ZT3,N);
        end
    end


    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A1)*FV(A2)*FV(A4)*valne(1,1)+...
                          FV(A1)*FV(A2)*DFV(A2)*DFV(A4)*(1-FV(A2))*valne(1,2)+...
                          FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A4)*valne(2,1)+...
                          FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A4)*(1-FV(A2))*valne(2,2));
                end
            end
        end
    end
end

function S4=case8(S4,N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,DF,SDEL,A1,A2,A3,A4,MIN)
    if N<4
        return
    end
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
    MIN=10^-4;
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

    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    FV=[F(I1),F(I2),F(I3),F(I4)];
                    DFV=[DF(I1),DF(I2),DF(I3),DF(I4)];
                    
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      valeq*(FV(A1)*FV(A2)*FV(A3)*FV(A4)*valne(1,1,1)+...
                           FV(A1)*FV(A2)*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,1,2)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(1,2,1)+...
                           FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(1,2,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*FV(A4)*valne(2,1,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,1,2)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*FV(A4)*valne(2,2,1)+...
                           FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*DFV(A3)*DFV(A4)*(1-FV(A3))*valne(2,2,2));
                end
            end
        end
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



