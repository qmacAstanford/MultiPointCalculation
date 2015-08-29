function S4=s4wlcnew(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d)
%S4 is a four point correlation function
%For example:
%   S4(1,1,1,1)=SAAAA
%   S4(1,1,2,1)=SAABA

S4=zeros(2,2,2,2);
MIN=5e-4;
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

    for M=0:ORDL-1
        for N1=1:ORDEig
            for N2=1:ORDEig
                for N3=1:ORDEig
                    
    % Case 1: J1==J2==J3==J4
        S4=case1(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,F,SDEL,MIN);
        S4=case1(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,F,SDEL,MIN);

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
        
        S4=case2(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z1,Z1L,F,DF,SDEL,1,2,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z1,Z1L,F,DF,SDEL,1,2,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z1,Z1L,F,DF,SDEL,1,3,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z1,Z1L,F,DF,SDEL,1,3,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z1,Z1L,F,DF,SDEL,1,4,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z1,Z1L,F,DF,SDEL,1,4,MIN);
                
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
    
        S4=case2(S4,N,NM,N1,N2,N3,M,R2,R12,R3,GL2,GLM12,GL3,YLM221,YLM,Z2,Z2L,F,DF,SDEL,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R2,R12,R4,GL2,GLM12,GL4,YLM221,YLM,Z2,Z2L,F,DF,SDEL,MIN);
        
        %....
        
        S4=case2(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z2,Z2L,F,DF,SDEL,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z2,Z2L,F,DF,SDEL,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z2,Z2L,F,DF,SDEL,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z2,Z2L,F,DF,SDEL,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z3,Z3L,F,DF,SDEL,MIN);
        S4=case2(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z3,Z3L,F,DF,SDEL,MIN);


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

function S4=case1(S4,N,NM,R1,R12,R3,F,MIN)
% Case 1: J1==J2==J3==J4
    E1=R1;
    E12=R12;
    E3=R3;
    
    % on same monomer
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
       valeq=(1/2).*E1.^(-4).*((-2).*(3+E1.*NM)+exp(1).^(E1.*NM).*(6+E1.*NM.*(( ...
         -4)+E1.*NM)));
    elseif (abs(E1-E12)<MIN && abs(E12-E3)>=MIN)
        valeq=E1.^(-3).*(E1+(-1).*E3).^(-2).*E3.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
          NM)).*E3.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E3.^2.*NM+E1.^3.*((-1)+ ...
          exp(1).^(E3.*NM)+(-1).*E3.*NM)+(-1).*E1.*E3.^2.*((-3)+E3.*NM+exp( ...
          1).^(E1.*NM).*(3+E3.*NM)));
    elseif (abs(E1-E3)<MIN && abs(E1-E12)>=MIN)
        valeq=E1.^(-3).*(E1+(-1).*E12).^(-2).*E12.^(-2).*(2.*((-1)+exp(1).^(E1.* ...
          NM)).*E12.^3+(2+exp(1).^(E1.*NM)).*E1.^2.*E12.^2.*NM+E1.^3.*((-1)+ ...
          exp(1).^(E12.*NM)+(-1).*E12.*NM)+(-1).*E1.*E12.^2.*((-3)+E12.*NM+ ...
          exp(1).^(E1.*NM).*(3+E12.*NM)));
    elseif (abs(E12-E3)<MIN && abs(E1-E3)>=MIN)
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
    valeq=finite(valeq);
    
    S4(1,1,1,1)=S4(1,1,1,1)+2*F(1)*N*valeq;
    S4(2,2,2,2)=S4(2,2,2,2)+2*F(2)*N*valeq;
end

function S4=case2(S4,N,NM,R1,R12,R3,Z1,Z1L,F,DF,SDEL,A1,A2,MIN)
% Case 2: J1 <J2==J3==J4
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=Z1;
    ZE1L=Z1L;
    
    % on same monomer
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
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
    
    % on different monomers
    valne1=(-1).*((-1)+ZE1).^(-2).*ZE1.*(1+N.*((-1)+ZE1)+(-1).*ZE1.^N);
    valne2=(-1).*((-1)+ZE1L).^(-2).*ZE1L.*(1+N.*((-1)+ZE1L)+(-1).*ZE1L.^N);
    valeq=finite(valeq);
    valne1=finite(valne1);
    valne2=finite(valne2);
    
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
    E1=R1;
    E12=R12;
    E3=R3;
    ZE12=Z12;
    ZE12L=Z12L;

    % on same monomer
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
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
    
    % on different monomers
    valne1=(-1).*((-1)+ZE12).^(-2).*ZE12.*(1+N.*((-1)+ZE12)+(-1).*ZE12.^N);
    valne2=(-1).*((-1)+ZE12L).^(-2).*ZE12L.*(1+N.*((-1)+ZE12L)+(-1).*ZE12L.^N);
    valeq=finite(valeq);
    valne1=finite(valne1);
    valne2=finite(valne2);
    
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
    E1=R1;
    E12=R12;
    E3=R3;
    ZE3=Z3;
    ZE3L=Z3L;
    
    % on same monomer
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
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
    end
    
    % on different monomers
    valne1=(-1).*((-1)+ZE3).^(-2).*ZE3.*(1+N.*((-1)+ZE3)+(-1).*ZE3.^N);
    valne2=(-1).*((-1)+ZE3L).^(-2).*ZE3L.*(1+N.*((-1)+ZE3L)+(-1).*ZE3L.^N);
    valeq=finite(valeq);
    valne1=finite(valne1);
    valne2=finite(valne2);
    
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

function S4=case5(S4,N,NM,R1,R12,R3,Z1,Z1L,Z12,Z12L,F,DF,SDEL,A1,A2,A3,MIN)
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];

    % on same monomer
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

    % on different monomers
    MIN=5e-3;
    valne=zeros(2,2);
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT12=ZE12(I2);
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
        end
    end
    valeq=finite(valeq);
    valne=finite(valne);

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
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];

    % on same monomer
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
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
       valeq=(-1).*(exp(1).^E12).^((-1).*NM).*(1+(-1).*(exp(1).^E3).^((-1).*NM) ...
       ).*((exp(1).^E12).^NM+(-1).*(exp(1).^E3).^NM).*log(exp(1).^E1).^( ...
       -1).*(log(exp(1).^E1)+(-1).*log(exp(1).^E12)).^(-1).*log(exp(1) ...
       .^E12).^(-1).*(((-1)+(exp(1).^E12).^NM).*log(exp(1).^E1)+(-1).*(( ...
       -1)+(exp(1).^E1).^NM).*log(exp(1).^E12)).*(log(exp(1).^E12)+(-1).* ...
       log(exp(1).^E3)).^(-1).*log(exp(1).^E3).^(-1);
    end
    
    % on different monomers
    MIN=5e-3;
    valne=zeros(2,2);
    for I1=1:2
        for I2=1:2
            ZT12=ZE12(I1);
            ZT3=ZE3(I2);
            if abs(ZT12-ZT3)<MIN
                valne(I1,I2)=((-1)+ZT12).^(-3).*ZT12.*(2.*ZT12+(-1).*N.*ZT12+N.*ZT12.^2+(-1).* ...
                  N.*ZT12.^N+(-2).*ZT12.^(1+N)+N.*ZT12.^(1+N));
            else
                valne(I1,I2)=((-1)+ZT12).^(-2).*ZT12.*(ZT12+(-1).*ZT3).^(-1).*((-1)+ZT3).^(-2) ...
                  .*ZT3.*((-2).*ZT12+N.*ZT12+ZT12.^2+(-1).*N.*ZT12.^2+ZT12.^N+2.* ...
                  ZT3+(-1).*N.*ZT3+N.*ZT12.^2.*ZT3+(-2).*ZT12.^N.*ZT3+(-1).*ZT3.^2+ ...
                  N.*ZT3.^2+(-1).*N.*ZT12.*ZT3.^2+ZT12.^N.*ZT3.^2+(-1).*ZT3.^N+2.* ...
                  ZT12.*ZT3.^N+(-1).*ZT12.^2.*ZT3.^N);
            end
        end
    end
    valeq=finite(valeq);
    valne=finite(valne);

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

function S4=case7(S4,N,NM,R1,R12,R3,Z1,Z1L,Z3,Z3L,F,DF,SDEL,A1,A2,A4,MIN)
    E1=R1;
    E12=R12;
    E3=R3;            
    ZE1=[Z1,Z1L];
    ZE3=[Z3,Z3L];

    % on same monomer
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
    
    % on different monomers
    MIN=5e-3;
    valne=zeros(2,2);
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT3=ZE3(I2);
            if abs(ZT1-ZT3)<MIN
                valne(I1,I2)=((-1)+ZT1).^(-3).*ZT1.*(2.*ZT1+(-1).*N.*ZT1+N.*ZT1.^2+(-1).*N.* ...
                  ZT1.^N+(-2).*ZT1.^(1+N)+N.*ZT1.^(1+N));
            else
                valne(I1,I2)=((-1)+ZT1).^(-2).*ZT1.*(ZT1+(-1).*ZT3).^(-1).*((-1)+ZT3).^(-2).* ...
                  ZT3.*((-2).*ZT1+N.*ZT1+ZT1.^2+(-1).*N.*ZT1.^2+ZT1.^N+2.*ZT3+(-1).* ...
                  N.*ZT3+N.*ZT1.^2.*ZT3+(-2).*ZT1.^N.*ZT3+(-1).*ZT3.^2+N.*ZT3.^2+( ...
                  -1).*N.*ZT1.*ZT3.^2+ZT1.^N.*ZT3.^2+(-1).*ZT3.^N+2.*ZT1.*ZT3.^N+( ...
                  -1).*ZT1.^2.*ZT3.^N);
            end
        end
    end
    valeq=finite(valeq);
    valne=finite(valne);

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
    E1=R1;
    E12=R12;
    E3=R3;
    ZE1=[Z1,Z1L];
    ZE12=[Z12,Z12L];
    ZE3=[Z3,Z3L];
    
    % on same monomer
    if (abs(E1-E12)<MIN && abs(E12-E3)<MIN)
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
       valeq=exp(1).^((-1).*(E1+E12+E3).*NM).*((-1)+exp(1).^(E1.*NM)).*(exp(1) ...
          .^(E1.*NM)+(-1).*exp(1).^(E12.*NM)).*(exp(1).^(E12.*NM)+(-1).*exp( ...
          1).^(E3.*NM)).*((-1)+exp(1).^(E3.*NM)).*E1.^(-1).*(E1+(-1).*E12) ...
          .^(-1).*(E12+(-1).*E3).^(-1).*E3.^(-1);
    end
    
    % on different monomers
    MIN=4e-3;
%     MIN=abs(exp(NM)*exp(R1)*(1-exp(MIN)));
%     MIN=1e-3;
    valne=zeros(2,2,2);
    for I1=1:2
        for I2=1:2
            for I3=1:2
                ZT1=ZE1(I1);
                ZT12=ZE12(I2);
                ZT3=ZE3(I3);

                % subcase 1
                if (abs(ZT1-1)<MIN && abs(ZT12-1)<MIN && abs(ZT3-1)<MIN)
                        valne(I1,I2,I3)=(1/24).*((-6).*N+11.*N.^2+(-6).*N.^3+N.^4);
                % subcase 2
                elseif (abs(ZT1-1)<MIN && abs(ZT12-1)<MIN && abs(ZT3-1)>=MIN)
                    valne(I1,I2,I3)=(-1/6).*((-1)+ZT3).^(-4).*ZT3.*(6+(-11).*N+6.*N.^2+(-1).*N.^3+18.* ...
                           N.*ZT3+(-15).*N.^2.*ZT3+3.*N.^3.*ZT3+(-9).*N.*ZT3.^2+12.*N.^2.* ...
                           ZT3.^2+(-3).*N.^3.*ZT3.^2+2.*N.*ZT3.^3+(-3).*N.^2.*ZT3.^3+N.^3.* ...
                           ZT3.^3+(-6).*ZT3.^N);
                % subcase 3
                elseif (abs(ZT1-1)>=MIN && abs(ZT12-1)<MIN && abs(ZT3-1)<MIN)
                    valne(I1,I2,I3)=(-1/6).*((-1)+ZT1).^(-4).*ZT1.*(12.*ZT1+(-4).*N.*ZT1+(-3).*N.^2.* ...
                          ZT1+N.^3.*ZT1+(-6).*ZT1.^2+2.*N.*ZT1.^2+6.*N.^2.*ZT1.^2+(-2).* ...
                         N.^3.*ZT1.^2+2.*N.*ZT1.^3+(-3).*N.^2.*ZT1.^3+N.^3.*ZT1.^3+(-6).* ...
                          ZT1.^N);
                % subcase 4
                elseif (abs(ZT1-1)<MIN && abs(ZT12-1)>=MIN && abs(ZT3-1)<MIN)
                    valne(I1,I2,I3)=(-1/6).*((-1)+ZT12).^(-4).*ZT12.*(6+(-11).*N+6.*N.^2+(-1).*N.^3+ ...
                   18.*N.*ZT12+(-15).*N.^2.*ZT12+3.*N.^3.*ZT12+(-9).*N.*ZT12.^2+12.* ...
                   N.^2.*ZT12.^2+(-3).*N.^3.*ZT12.^2+2.*N.*ZT12.^3+(-3).*N.^2.* ...
                   ZT12.^3+N.^3.*ZT12.^3+(-6).*ZT12.^N);
                % subcase 5
                elseif (abs(ZT1-1)<MIN && abs(ZT12-1)>=MIN && abs(ZT3-1)>=MIN)

                    if (abs(ZT12-ZT3)<MIN)
                        valne(I1,I2,I3)=(1/2).*((-1)+ZT12).^(-4).*ZT12.*(6.*ZT12+(-5).*N.*ZT12+N.^2.*ZT12+ ...
                          6.*N.*ZT12.^2+(-2).*N.^2.*ZT12.^2+(-1).*N.*ZT12.^3+N.^2.*ZT12.^3+( ...
                          -2).*N.*ZT12.^N+(-6).*ZT12.^(1+N)+2.*N.*ZT12.^(1+N));
                    else
                        valne(I1,I2,I3)=(1/2).*((-1)+ZT12).^(-3).*(ZT12+(-1).*ZT3).^(-1).*((-1)+ZT3).^(-3) ...
                          .*(6.*ZT12.^2.*ZT3+(-5).*N.*ZT12.^2.*ZT3+N.^2.*ZT12.^2.*ZT3+(-6).* ...
                          ZT12.^3.*ZT3+8.*N.*ZT12.^3.*ZT3+(-2).*N.^2.*ZT12.^3.*ZT3+2.* ...
                          ZT12.^4.*ZT3+(-3).*N.*ZT12.^4.*ZT3+N.^2.*ZT12.^4.*ZT3+(-2).* ...
                          ZT12.^(1+N).*ZT3+(-6).*ZT12.*ZT3.^2+5.*N.*ZT12.*ZT3.^2+(-1).* ...
                          N.^2.*ZT12.*ZT3.^2+(-9).*N.*ZT12.^3.*ZT3.^2+3.*N.^2.*ZT12.^3.* ...
                          ZT3.^2+4.*N.*ZT12.^4.*ZT3.^2+(-2).*N.^2.*ZT12.^4.*ZT3.^2+6.* ...
                          ZT12.^(1+N).*ZT3.^2+6.*ZT12.*ZT3.^3+(-8).*N.*ZT12.*ZT3.^3+2.* ...
                          N.^2.*ZT12.*ZT3.^3+9.*N.*ZT12.^2.*ZT3.^3+(-3).*N.^2.*ZT12.^2.* ...
                          ZT3.^3+(-1).*N.*ZT12.^4.*ZT3.^3+N.^2.*ZT12.^4.*ZT3.^3+(-6).* ...
                          ZT12.^(1+N).*ZT3.^3+(-2).*ZT12.*ZT3.^4+3.*N.*ZT12.*ZT3.^4+(-1).* ...
                          N.^2.*ZT12.*ZT3.^4+(-4).*N.*ZT12.^2.*ZT3.^4+2.*N.^2.*ZT12.^2.* ...
                          ZT3.^4+N.*ZT12.^3.*ZT3.^4+(-1).*N.^2.*ZT12.^3.*ZT3.^4+2.*ZT12.^(1+ ...
                          N).*ZT3.^4+2.*ZT12.*ZT3.^(1+N)+(-6).*ZT12.^2.*ZT3.^(1+N)+6.* ...
                          ZT12.^3.*ZT3.^(1+N)+(-2).*ZT12.^4.*ZT3.^(1+N));
                    end
                % subcase 6
                elseif (abs(ZT1-1)>=MIN && abs(ZT12-1)<MIN && abs(ZT3-1)>=MIN)
                    
                    if (abs(ZT1-ZT3)<MIN)
                        valne(I1,I2,I3)=(1/2).*((-1)+ZT1).^(-4).*ZT1.*(6.*ZT1+(-5).*N.*ZT1+N.^2.*ZT1+6.* ...
                          N.*ZT1.^2+(-2).*N.^2.*ZT1.^2+(-1).*N.*ZT1.^3+N.^2.*ZT1.^3+(-2).* ...
                          N.*ZT1.^N+(-6).*ZT1.^(1+N)+2.*N.*ZT1.^(1+N));
                    else
                        valne(I1,I2,I3)=(1/2).*((-1)+ZT1).^(-3).*(ZT1+(-1).*ZT3).^(-1).*((-1)+ZT3).^(-3).* ...
                          (6.*ZT1.^2.*ZT3+(-5).*N.*ZT1.^2.*ZT3+N.^2.*ZT1.^2.*ZT3+(-6).* ...
                          ZT1.^3.*ZT3+8.*N.*ZT1.^3.*ZT3+(-2).*N.^2.*ZT1.^3.*ZT3+2.*ZT1.^4.* ...
                          ZT3+(-3).*N.*ZT1.^4.*ZT3+N.^2.*ZT1.^4.*ZT3+(-2).*ZT1.^(1+N).*ZT3+( ...
                          -6).*ZT1.*ZT3.^2+5.*N.*ZT1.*ZT3.^2+(-1).*N.^2.*ZT1.*ZT3.^2+(-9).* ...
                          N.*ZT1.^3.*ZT3.^2+3.*N.^2.*ZT1.^3.*ZT3.^2+4.*N.*ZT1.^4.*ZT3.^2+( ...
                          -2).*N.^2.*ZT1.^4.*ZT3.^2+6.*ZT1.^(1+N).*ZT3.^2+6.*ZT1.*ZT3.^3+( ...
                          -8).*N.*ZT1.*ZT3.^3+2.*N.^2.*ZT1.*ZT3.^3+9.*N.*ZT1.^2.*ZT3.^3+(-3) ...
                          .*N.^2.*ZT1.^2.*ZT3.^3+(-1).*N.*ZT1.^4.*ZT3.^3+N.^2.*ZT1.^4.* ...
                          ZT3.^3+(-6).*ZT1.^(1+N).*ZT3.^3+(-2).*ZT1.*ZT3.^4+3.*N.*ZT1.* ...
                          ZT3.^4+(-1).*N.^2.*ZT1.*ZT3.^4+(-4).*N.*ZT1.^2.*ZT3.^4+2.*N.^2.* ...
                          ZT1.^2.*ZT3.^4+N.*ZT1.^3.*ZT3.^4+(-1).*N.^2.*ZT1.^3.*ZT3.^4+2.* ...
                          ZT1.^(1+N).*ZT3.^4+2.*ZT1.*ZT3.^(1+N)+(-6).*ZT1.^2.*ZT3.^(1+N)+6.* ...
                          ZT1.^3.*ZT3.^(1+N)+(-2).*ZT1.^4.*ZT3.^(1+N));
                    end
                % subcase 7
                elseif (abs(ZT1-1)>=MIN && abs(ZT12-1)>=MIN && abs(ZT3-1)<MIN)
                    
                    if (abs(ZT1-ZT12)<MIN)
                        valne(I1,I2,I3)=(1/2).*((-1)+ZT1).^(-4).*ZT1.*(6.*ZT1+(-5).*N.*ZT1+N.^2.*ZT1+6.* ...
                           N.*ZT1.^2+(-2).*N.^2.*ZT1.^2+(-1).*N.*ZT1.^3+N.^2.*ZT1.^3+(-2).* ...
                           N.*ZT1.^N+(-6).*ZT1.^(1+N)+2.*N.*ZT1.^(1+N));
                    else
                        valne(I1,I2,I3)=(1/2).*((-1)+ZT1).^(-3).*(ZT1+(-1).*ZT12).^(-1).*((-1)+ZT12).^(-3) ...
                          .*(6.*ZT1.^2.*ZT12+(-5).*N.*ZT1.^2.*ZT12+N.^2.*ZT1.^2.*ZT12+(-6).* ...
                          ZT1.^3.*ZT12+8.*N.*ZT1.^3.*ZT12+(-2).*N.^2.*ZT1.^3.*ZT12+2.* ...
                          ZT1.^4.*ZT12+(-3).*N.*ZT1.^4.*ZT12+N.^2.*ZT1.^4.*ZT12+(-2).*ZT1.^( ...
                          1+N).*ZT12+(-6).*ZT1.*ZT12.^2+5.*N.*ZT1.*ZT12.^2+(-1).*N.^2.*ZT1.* ...
                          ZT12.^2+(-9).*N.*ZT1.^3.*ZT12.^2+3.*N.^2.*ZT1.^3.*ZT12.^2+4.*N.* ...
                          ZT1.^4.*ZT12.^2+(-2).*N.^2.*ZT1.^4.*ZT12.^2+6.*ZT1.^(1+N).* ...
                          ZT12.^2+6.*ZT1.*ZT12.^3+(-8).*N.*ZT1.*ZT12.^3+2.*N.^2.*ZT1.* ...
                          ZT12.^3+9.*N.*ZT1.^2.*ZT12.^3+(-3).*N.^2.*ZT1.^2.*ZT12.^3+(-1).* ...
                          N.*ZT1.^4.*ZT12.^3+N.^2.*ZT1.^4.*ZT12.^3+(-6).*ZT1.^(1+N).* ...
                          ZT12.^3+(-2).*ZT1.*ZT12.^4+3.*N.*ZT1.*ZT12.^4+(-1).*N.^2.*ZT1.* ...
                          ZT12.^4+(-4).*N.*ZT1.^2.*ZT12.^4+2.*N.^2.*ZT1.^2.*ZT12.^4+N.* ...
                          ZT1.^3.*ZT12.^4+(-1).*N.^2.*ZT1.^3.*ZT12.^4+2.*ZT1.^(1+N).* ...
                          ZT12.^4+2.*ZT1.*ZT12.^(1+N)+(-6).*ZT1.^2.*ZT12.^(1+N)+6.*ZT1.^3.* ...
                          ZT12.^(1+N)+(-2).*ZT1.^4.*ZT12.^(1+N));
                    end
                %subcase 8
                else
                    
                    if (abs(ZT1-ZT12)<MIN && abs(ZT12-ZT3)<MIN)
                        valne(I1,I2,I3)=(1/2).*((-1)+ZT1).^(-4).*(2.*ZT1.^3.*((-3)+N+(-1).*N.*ZT1)+ZT1.^( ...
                      1+N).*(((-1)+N).*N+(-2).*((-3)+N).*N.*ZT1+((-3)+N).*((-2)+N).* ...
                     ZT1.^2));
                    elseif (abs(ZT1-ZT12)<MIN && abs(ZT1-ZT3)>=MIN)
                        valne(I1,I2,I3)=((-1)+ZT1).^(-3).*(ZT1+(-1).*ZT3).^(-2).*((-1)+ZT3).^(-2).*(3.* ...
                       ZT1.^4.*ZT3+(-1).*N.*ZT1.^4.*ZT3+(-1).*ZT1.^5.*ZT3+N.*ZT1.^5.*ZT3+ ...
                       ZT1.^(2+N).*ZT3+(-1).*N.*ZT1.^(2+N).*ZT3+(-3).*ZT1.^(3+N).*ZT3+N.* ...
                       ZT1.^(3+N).*ZT3+(-6).*ZT1.^3.*ZT3.^2+2.*N.*ZT1.^3.*ZT3.^2+(-1).* ...
                       N.*ZT1.^4.*ZT3.^2+(-1).*N.*ZT1.^5.*ZT3.^2+N.*ZT1.^(1+N).*ZT3.^2+ ...
                       N.*ZT1.^(2+N).*ZT3.^2+6.*ZT1.^(3+N).*ZT3.^2+(-2).*N.*ZT1.^(3+N).* ...
                       ZT3.^2+3.*ZT1.^2.*ZT3.^3+(-1).*N.*ZT1.^2.*ZT3.^3+3.*ZT1.^3.* ...
                       ZT3.^3+(-1).*N.*ZT1.^3.*ZT3.^3+2.*N.*ZT1.^4.*ZT3.^3+(-2).*N.* ...
                       ZT1.^(1+N).*ZT3.^3+(-3).*ZT1.^(2+N).*ZT3.^3+N.*ZT1.^(2+N).*ZT3.^3+ ...
                       (-3).*ZT1.^(3+N).*ZT3.^3+N.*ZT1.^(3+N).*ZT3.^3+(-2).*ZT1.^2.* ...
                       ZT3.^4+N.*ZT1.^2.*ZT3.^4+(-1).*N.*ZT1.^3.*ZT3.^4+N.*ZT1.^(1+N).* ...
                       ZT3.^4+2.*ZT1.^(2+N).*ZT3.^4+(-1).*N.*ZT1.^(2+N).*ZT3.^4+(-1).* ...
                       ZT1.^2.*ZT3.^(1+N)+3.*ZT1.^3.*ZT3.^(1+N)+(-3).*ZT1.^4.*ZT3.^(1+N)+ ...
                       ZT1.^5.*ZT3.^(1+N));
                    elseif (abs(ZT1-ZT3)<MIN && abs(ZT1-ZT12)>=MIN)
                        valne(I1,I2,I3)=((-1)+ZT1).^(-3).*(ZT1+(-1).*ZT12).^(-1).*((-1)+ZT12).^(-2).*((-1) ...
                      .*ZT1+ZT12).^(-1).*((-3).*ZT1.^4.*ZT12+N.*ZT1.^4.*ZT12+ZT1.^5.* ...
                      ZT12+(-1).*N.*ZT1.^5.*ZT12+(-1).*ZT1.^(2+N).*ZT12+N.*ZT1.^(2+N).* ...
                      ZT12+3.*ZT1.^(3+N).*ZT12+(-1).*N.*ZT1.^(3+N).*ZT12+6.*ZT1.^3.* ...
                      ZT12.^2+(-2).*N.*ZT1.^3.*ZT12.^2+N.*ZT1.^4.*ZT12.^2+N.*ZT1.^5.* ...
                      ZT12.^2+(-1).*N.*ZT1.^(1+N).*ZT12.^2+(-1).*N.*ZT1.^(2+N).*ZT12.^2+ ...
                      (-6).*ZT1.^(3+N).*ZT12.^2+2.*N.*ZT1.^(3+N).*ZT12.^2+(-3).*ZT1.^2.* ...
                      ZT12.^3+N.*ZT1.^2.*ZT12.^3+(-3).*ZT1.^3.*ZT12.^3+N.*ZT1.^3.* ...
                      ZT12.^3+(-2).*N.*ZT1.^4.*ZT12.^3+2.*N.*ZT1.^(1+N).*ZT12.^3+3.* ...
                      ZT1.^(2+N).*ZT12.^3+(-1).*N.*ZT1.^(2+N).*ZT12.^3+3.*ZT1.^(3+N).* ...
                      ZT12.^3+(-1).*N.*ZT1.^(3+N).*ZT12.^3+2.*ZT1.^2.*ZT12.^4+(-1).*N.* ...
                      ZT1.^2.*ZT12.^4+N.*ZT1.^3.*ZT12.^4+(-1).*N.*ZT1.^(1+N).*ZT12.^4+( ...
                      -2).*ZT1.^(2+N).*ZT12.^4+N.*ZT1.^(2+N).*ZT12.^4+ZT1.^2.*ZT12.^(1+ ...
                      N)+(-3).*ZT1.^3.*ZT12.^(1+N)+3.*ZT1.^4.*ZT12.^(1+N)+(-1).*ZT1.^5.* ...
                      ZT12.^(1+N));
                    elseif (abs(ZT12-ZT3)<MIN &&  abs(ZT1-ZT3)>=MIN)
                        valne(I1,I2,I3)=(-1).*((-1)+ZT1).^(-2).*(ZT1+(-1).*ZT12).^(-2).*((-1)+ZT12).^(-3) ...
                      .*((-3).*ZT1.^3.*ZT12.^2+N.*ZT1.^3.*ZT12.^2+2.*ZT1.^4.*ZT12.^2+( ...
                      -1).*N.*ZT1.^4.*ZT12.^2+ZT1.^(1+N).*ZT12.^2+6.*ZT1.^2.*ZT12.^3+( ...
                      -2).*N.*ZT1.^2.*ZT12.^3+(-3).*ZT1.^3.*ZT12.^3+N.*ZT1.^3.*ZT12.^3+ ...
                      N.*ZT1.^4.*ZT12.^3+(-3).*ZT1.^(1+N).*ZT12.^3+(-3).*ZT1.*ZT12.^4+ ...
                      N.*ZT1.*ZT12.^4+N.*ZT1.^2.*ZT12.^4+(-2).*N.*ZT1.^3.*ZT12.^4+3.* ...
                      ZT1.^(1+N).*ZT12.^4+ZT1.*ZT12.^5+(-1).*N.*ZT1.*ZT12.^5+N.*ZT1.^2.* ...
                      ZT12.^5+(-1).*ZT1.^(1+N).*ZT12.^5+(-1).*N.*ZT1.^2.*ZT12.^(1+N)+2.* ...
                      N.*ZT1.^3.*ZT12.^(1+N)+(-1).*N.*ZT1.^4.*ZT12.^(1+N)+(-1).*ZT1.* ...
                      ZT12.^(2+N)+N.*ZT1.*ZT12.^(2+N)+(-1).*N.*ZT1.^2.*ZT12.^(2+N)+3.* ...
                      ZT1.^3.*ZT12.^(2+N)+(-1).*N.*ZT1.^3.*ZT12.^(2+N)+(-2).*ZT1.^4.* ...
                      ZT12.^(2+N)+N.*ZT1.^4.*ZT12.^(2+N)+3.*ZT1.*ZT12.^(3+N)+(-1).*N.* ...
                      ZT1.*ZT12.^(3+N)+(-6).*ZT1.^2.*ZT12.^(3+N)+2.*N.*ZT1.^2.*ZT12.^(3+ ...
                      N)+3.*ZT1.^3.*ZT12.^(3+N)+(-1).*N.*ZT1.^3.*ZT12.^(3+N));
                    else
                        valne(I1,I2,I3)=((-1)+ZT1).^(-2).*(ZT1+(-1).*ZT12).^(-1).*((-1)+ZT12).^(-2).*(ZT1+ ...
                       (-1).*ZT3).^(-1).*((-1)+ZT3).^(-2).*((-1).*ZT12+ZT3).^(-1).*(3.* ...
                       ZT1.^3.*ZT12.^2.*ZT3+(-1).*N.*ZT1.^3.*ZT12.^2.*ZT3+(-2).*ZT1.^4.* ...
                       ZT12.^2.*ZT3+N.*ZT1.^4.*ZT12.^2.*ZT3+(-1).*ZT1.^(1+N).*ZT12.^2.* ...
                       ZT3+(-3).*ZT1.^2.*ZT12.^3.*ZT3+N.*ZT1.^2.*ZT12.^3.*ZT3+ZT1.^4.* ...
                       ZT12.^3.*ZT3+(-1).*N.*ZT1.^4.*ZT12.^3.*ZT3+2.*ZT1.^(1+N).* ...
                       ZT12.^3.*ZT3+2.*ZT1.^2.*ZT12.^4.*ZT3+(-1).*N.*ZT1.^2.*ZT12.^4.* ...
                       ZT3+(-1).*ZT1.^3.*ZT12.^4.*ZT3+N.*ZT1.^3.*ZT12.^4.*ZT3+(-1).* ...
                       ZT1.^(1+N).*ZT12.^4.*ZT3+ZT1.^2.*ZT12.^(1+N).*ZT3+(-2).*ZT1.^3.* ...
                       ZT12.^(1+N).*ZT3+ZT1.^4.*ZT12.^(1+N).*ZT3+(-3).*ZT1.^3.*ZT12.* ...
                       ZT3.^2+N.*ZT1.^3.*ZT12.*ZT3.^2+2.*ZT1.^4.*ZT12.*ZT3.^2+(-1).*N.* ...
                       ZT1.^4.*ZT12.*ZT3.^2+ZT1.^(1+N).*ZT12.*ZT3.^2+3.*ZT1.*ZT12.^3.* ...
                       ZT3.^2+(-1).*N.*ZT1.*ZT12.^3.*ZT3.^2+N.*ZT1.^4.*ZT12.^3.*ZT3.^2+( ...
                       -3).*ZT1.^(1+N).*ZT12.^3.*ZT3.^2+(-2).*ZT1.*ZT12.^4.*ZT3.^2+N.* ...
                       ZT1.*ZT12.^4.*ZT3.^2+(-1).*N.*ZT1.^3.*ZT12.^4.*ZT3.^2+2.*ZT1.^(1+ ...
                       N).*ZT12.^4.*ZT3.^2+(-1).*ZT1.*ZT12.^(1+N).*ZT3.^2+3.*ZT1.^3.* ...
                       ZT12.^(1+N).*ZT3.^2+(-2).*ZT1.^4.*ZT12.^(1+N).*ZT3.^2+3.*ZT1.^2.* ...
                       ZT12.*ZT3.^3+(-1).*N.*ZT1.^2.*ZT12.*ZT3.^3+(-1).*ZT1.^4.*ZT12.* ...
                       ZT3.^3+N.*ZT1.^4.*ZT12.*ZT3.^3+(-2).*ZT1.^(1+N).*ZT12.*ZT3.^3+(-3) ...
                       .*ZT1.*ZT12.^2.*ZT3.^3+N.*ZT1.*ZT12.^2.*ZT3.^3+(-1).*N.*ZT1.^4.* ...
                       ZT12.^2.*ZT3.^3+3.*ZT1.^(1+N).*ZT12.^2.*ZT3.^3+ZT1.*ZT12.^4.* ...
                       ZT3.^3+(-1).*N.*ZT1.*ZT12.^4.*ZT3.^3+N.*ZT1.^2.*ZT12.^4.*ZT3.^3+( ...
                       -1).*ZT1.^(1+N).*ZT12.^4.*ZT3.^3+2.*ZT1.*ZT12.^(1+N).*ZT3.^3+(-3) ...
                       .*ZT1.^2.*ZT12.^(1+N).*ZT3.^3+ZT1.^4.*ZT12.^(1+N).*ZT3.^3+(-2).* ...
                       ZT1.^2.*ZT12.*ZT3.^4+N.*ZT1.^2.*ZT12.*ZT3.^4+ZT1.^3.*ZT12.*ZT3.^4+ ...
                       (-1).*N.*ZT1.^3.*ZT12.*ZT3.^4+ZT1.^(1+N).*ZT12.*ZT3.^4+2.*ZT1.* ...
                       ZT12.^2.*ZT3.^4+(-1).*N.*ZT1.*ZT12.^2.*ZT3.^4+N.*ZT1.^3.*ZT12.^2.* ...
                       ZT3.^4+(-2).*ZT1.^(1+N).*ZT12.^2.*ZT3.^4+(-1).*ZT1.*ZT12.^3.* ...
                       ZT3.^4+N.*ZT1.*ZT12.^3.*ZT3.^4+(-1).*N.*ZT1.^2.*ZT12.^3.*ZT3.^4+ ...
                       ZT1.^(1+N).*ZT12.^3.*ZT3.^4+(-1).*ZT1.*ZT12.^(1+N).*ZT3.^4+2.* ...
                       ZT1.^2.*ZT12.^(1+N).*ZT3.^4+(-1).*ZT1.^3.*ZT12.^(1+N).*ZT3.^4+(-1) ...
                       .*ZT1.^2.*ZT12.*ZT3.^(1+N)+2.*ZT1.^3.*ZT12.*ZT3.^(1+N)+(-1).* ...
                       ZT1.^4.*ZT12.*ZT3.^(1+N)+ZT1.*ZT12.^2.*ZT3.^(1+N)+(-3).*ZT1.^3.* ...
                       ZT12.^2.*ZT3.^(1+N)+2.*ZT1.^4.*ZT12.^2.*ZT3.^(1+N)+(-2).*ZT1.* ...
                       ZT12.^3.*ZT3.^(1+N)+3.*ZT1.^2.*ZT12.^3.*ZT3.^(1+N)+(-1).*ZT1.^4.* ...
                       ZT12.^3.*ZT3.^(1+N)+ZT1.*ZT12.^4.*ZT3.^(1+N)+(-2).*ZT1.^2.* ...
                       ZT12.^4.*ZT3.^(1+N)+ZT1.^3.*ZT12.^4.*ZT3.^(1+N));
                    end
                end
                
%                 if (valne(I1,I2,I3)>1e5)
%                     token
%                 end
            end
        end
    end
%     valne=ones(2,2,2);
    valeq=finite(valeq);
    valne=finite(valne);
    
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

function xf=finite(x)
    MIN=-1e5;
    MAX=1e5;
    xf=x;
    xf(xf>MAX | xf<MIN)=0;
    xf(isnan(xf))=0;
end