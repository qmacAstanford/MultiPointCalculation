function S4=s4wlc(N,NM,LAM,FA,Q1,Q2,Q3,Q4,ORDEig,ORDL,NumLayer)

S4=zeros(2,2,2,2);
MIN=5e-4;

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
    
    QMIN=3e-1/N/NM;
    Q1MAG=max([Q1MAG,QMIN]);
    Q2MAG=max([Q2MAG,QMIN]);
    Q3MAG=max([Q3MAG,QMIN]);
    Q4MAG=max([Q4MAG,QMIN]);
    
    Q12MAG=sqrt(sum(power(Q1+Q2,2)));
    Q13MAG=sqrt(sum(power(Q1+Q3,2)));
    Q14MAG=sqrt(sum(power(Q1+Q4,2)));
    Q23MAG=sqrt(sum(power(Q2+Q3,2)));
    Q24MAG=sqrt(sum(power(Q2+Q4,2)));
    Q34MAG=sqrt(sum(power(Q3+Q4,2)));    

    Q12MAG=max([Q12MAG,QMIN]);
    Q13MAG=max([Q13MAG,QMIN]);
    Q14MAG=max([Q14MAG,QMIN]);
    Q23MAG=max([Q23MAG,QMIN]);
    Q24MAG=max([Q24MAG,QMIN]);
    Q34MAG=max([Q34MAG,QMIN]);
    
    EQ1=Q1/Q1MAG;
    EQ2=Q2/Q2MAG;
    EQ3=Q3/Q3MAG;
    EQ4=Q4/Q4MAG;
    
    EQ12=(Q1+Q2)/Q12MAG;
    EQ14=(Q1+Q4)/Q14MAG;
    EQ23=(Q2+Q3)/Q23MAG;
    EQ34=(Q3+Q4)/Q34MAG;
    
    [RHO112,PHI112]=rotvec(EQ1,EQ12);
    [RHO114,PHI114]=rotvec(EQ1,EQ14);
    [RHO221,PHI221]=rotvec(EQ2,EQ12);
    [RHO223,PHI223]=rotvec(EQ2,EQ23);
    [RHO332,PHI332]=rotvec(EQ3,EQ23);
    [RHO334,PHI334]=rotvec(EQ3,EQ34);
    [RHO441,PHI441]=rotvec(EQ4,EQ14);
    [RHO443,PHI443]=rotvec(EQ4,EQ34);
    
    YLM112=zeros(ORDL,ORDL);
    YLM113=zeros(ORDL,ORDL);
    YLM114=zeros(ORDL,ORDL);
    YLM221=zeros(ORDL,ORDL);
    YLM223=zeros(ORDL,ORDL);
    YLM224=zeros(ORDL,ORDL);
    YLM331=zeros(ORDL,ORDL);
    YLM332=zeros(ORDL,ORDL); 
    YLM334=zeros(ORDL,ORDL);
    YLM441=zeros(ORDL,ORDL);
    YLM442=zeros(ORDL,ORDL);
    YLM443=zeros(ORDL,ORDL);
    for L=0:(ORDL-1)
        PLM112=legendre(L,RHO112);
        PLM114=legendre(L,RHO114);
        PLM221=legendre(L,RHO221);
        PLM223=legendre(L,RHO223);
        PLM332=legendre(L,RHO332);
        PLM334=legendre(L,RHO334);
        PLM441=legendre(L,RHO441);
        PLM443=legendre(L,RHO443);
        for M=0:L
            YLM112(L+1,M+1)=exp(1i*M*PHI112)*sqrt(factorial(L-M)/factorial(L+M))*PLM112(M+1);
            YLM114(L+1,M+1)=exp(1i*M*PHI114)*sqrt(factorial(L-M)/factorial(L+M))*PLM114(M+1);
            YLM221(L+1,M+1)=exp(1i*M*PHI221)*sqrt(factorial(L-M)/factorial(L+M))*PLM221(M+1);
            YLM223(L+1,M+1)=exp(1i*M*PHI223)*sqrt(factorial(L-M)/factorial(L+M))*PLM223(M+1);
            YLM332(L+1,M+1)=exp(1i*M*PHI332)*sqrt(factorial(L-M)/factorial(L+M))*PLM332(M+1);
            YLM334(L+1,M+1)=exp(1i*M*PHI334)*sqrt(factorial(L-M)/factorial(L+M))*PLM334(M+1);
            YLM441(L+1,M+1)=exp(1i*M*PHI441)*sqrt(factorial(L-M)/factorial(L+M))*PLM441(M+1);
            YLM443(L+1,M+1)=exp(1i*M*PHI443)*sqrt(factorial(L-M)/factorial(L+M))*PLM443(M+1);
        end
    end
    
    R1=MatRoots(Q1MAG,3,ORDEig);
    R2=MatRoots(Q2MAG,3,ORDEig);
    R3=MatRoots(Q3MAG,3,ORDEig);
    R4=MatRoots(Q4MAG,3,ORDEig);
    
    R12=MatRootsLM(Q12MAG,ORDEig);
    R13=MatRootsLM(Q13MAG,ORDEig);
    R14=MatRootsLM(Q14MAG,ORDEig);
    R23=MatRootsLM(Q23MAG,ORDEig);
    R24=MatRootsLM(Q24MAG,ORDEig);
    R34=MatRootsLM(Q34MAG,ORDEig);
    
    GL1=gl(Q1MAG,R1,ORDEig,ORDL,NumLayer);
    GL2=gl(Q2MAG,R2,ORDEig,ORDL,NumLayer);
    GL3=gl(Q3MAG,R3,ORDEig,ORDL,NumLayer);
    GL4=gl(Q3MAG,R4,ORDEig,ORDL,NumLayer);
    
    GLM12=glm(Q12MAG,R12,ORDEig,ORDL,NumLayer);
    GLM13=glm(Q13MAG,R13,ORDEig,ORDL,NumLayer);
    GLM14=glm(Q14MAG,R14,ORDEig,ORDL,NumLayer);
    GLM23=glm(Q23MAG,R23,ORDEig,ORDL,NumLayer);
    GLM24=glm(Q24MAG,R24,ORDEig,ORDL,NumLayer);
    GLM34=glm(Q34MAG,R34,ORDEig,ORDL,NumLayer);

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

    for M=0:ORDL-1
        for N1=1:ORDEig
            for N2=1:ORDEig
                for N3=1:ORDEig

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

                % Case 1: J1==J2==J3==J4
                SDEL=zeros(2,2,2,2);
                SDEL(1,1,1,1)=1;
                SDEL(2,2,2,2)=1;

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

                % Case 2: J1 <J2==J3==J4
                SDEL=zeros(2,2,2,2);
                SDEL(1,1,1,1)=1;
                SDEL(2,2,2,2)=1;
                SDEL(2,1,1,1)=1;
                SDEL(1,2,2,2)=1;

                S4=case2(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case2(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z3,Z3L,F,DF,SDEL,MIN);

                % Case 3: J1==J2 <J3==J4
                SDEL=zeros(2,2,2,2);
                SDEL(1,1,1,1)=1;
                SDEL(2,2,2,2)=1;
                SDEL(1,1,2,2)=1;
                SDEL(2,2,1,1)=1;

                S4=case3(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z12,Z12L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z12,Z12L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z13,Z13L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z13,Z13L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z14,Z14L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z14,Z14L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z23,Z23L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z23,Z23L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z24,Z24L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z24,Z24L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z34,Z34L,F,DF,SDEL,MIN);
                S4=case3(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z34,Z34L,F,DF,SDEL,MIN);
                
                % Case 4: J1==J2==J3 <J4
                SDEL=zeros(2,2,2,2);
                SDEL(1,1,1,1)=1;
                SDEL(2,2,2,2)=1;
                SDEL(1,1,1,2)=1;
                SDEL(2,2,2,1)=1;

                S4=case4(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case4(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z1,Z1L,F,DF,SDEL,MIN);

                % Case 5: J1 <J2 <J3==J4
                SDEL=zeros(2,2,2,2);
                SDEL(1,1,1,1)=1;
                SDEL(2,2,2,2)=1;
                SDEL(2,1,1,1)=1;
                SDEL(1,2,1,1)=1;
                SDEL(2,2,1,1)=1;
                SDEL(1,2,2,2)=1;
                SDEL(2,1,2,2)=1;
                SDEL(1,1,2,2)=1;

                S4=case5(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z1,Z1L,Z12,Z12L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z1,Z1L,Z12,Z12L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z1,Z1L,Z13,Z13L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z1,Z1L,Z13,Z13L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z1,Z1L,Z14,Z14L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z1,Z1L,Z14,Z14L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z2,Z2L,Z23,Z23L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z2,Z2L,Z23,Z23L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z2,Z2L,Z24,Z24L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z2,Z2L,Z24,Z24L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z3,Z3L,Z34,Z34L,F,DF,SDEL,MIN);
                S4=case5(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z3,Z3L,Z34,Z34L,F,DF,SDEL,MIN);

                % Case 6: J1==J2 <J3 <J4
                SDEL=zeros(2,2,2,2);
                SDEL(1,1,1,1)=1;
                SDEL(2,2,2,2)=1;
                SDEL(1,1,1,2)=1;
                SDEL(1,1,2,1)=1;
                SDEL(1,1,2,2)=1;
                SDEL(2,2,2,1)=1;
                SDEL(2,2,1,2)=1;
                SDEL(2,2,1,1)=1;

                S4=case6(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z12,Z12L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z12,Z12L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z13,Z13L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z13,Z13L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z14,Z14L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z14,Z14L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z23,Z23L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z23,Z23L,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z24,Z24L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z24,Z24L,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z34,Z34L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case6(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z34,Z34L,Z1,Z1L,F,DF,SDEL,MIN);

                % Case 7: J1 <J2==J3 <J4
                SDEL=zeros(2,2,2,2);
                SDEL(1,1,1,1)=1;
                SDEL(2,2,2,2)=1;
                SDEL(1,1,1,2)=1;
                SDEL(2,1,1,1)=1;
                SDEL(2,1,1,2)=1;
                SDEL(2,2,2,1)=1;
                SDEL(1,2,2,2)=1;
                SDEL(1,2,2,1)=1;

                S4=case7(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z1,Z1L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z1,Z1L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z1,Z1L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z1,Z1L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z1,Z1L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z1,Z1L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z2,Z2L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z2,Z2L,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z2,Z2L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z2,Z2L,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z3,Z3L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case7(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z3,Z3L,Z1,Z1L,F,DF,SDEL,MIN);
                
                % Case 8: J1 <J2 <J3 <J4
                SDEL=ones(2,2,2,2);

                S4=case8(S4,N,NM,N1,N2,N3,M,R1,R12,R4,GL1,GLM12,GL4,YLM112,YLM443,Z1,Z1L,Z12,Z12L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM334,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R1,R13,R4,GL1,GLM13,GL4,YLM113,YLM442,Z1,Z1L,Z13,Z13L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R1,R13,R2,GL1,GLM13,GL2,YLM113,YLM224,Z1,Z1L,Z13,Z13L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R1,R14,R3,GL1,GLM14,GL3,YLM114,YLM332,Z1,Z1L,Z14,Z14L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R1,R14,R2,GL1,GLM14,GL2,YLM114,YLM223,Z1,Z1L,Z14,Z14L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R2,R23,R4,GL2,GLM23,GL4,YLM223,YLM441,Z2,Z2L,Z23,Z23L,Z4,Z4L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R2,R23,R1,GL2,GLM23,GL1,YLM223,YLM114,Z2,Z2L,Z23,Z23L,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R2,R24,R3,GL2,GLM24,GL3,YLM224,YLM331,Z2,Z2L,Z24,Z24L,Z3,Z3L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R2,R24,R1,GL2,GLM24,GL1,YLM224,YLM113,Z2,Z2L,Z24,Z24L,Z1,Z1L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R3,R34,R2,GL3,GLM34,GL2,YLM334,YLM221,Z3,Z3L,Z34,Z34L,Z2,Z2L,F,DF,SDEL,MIN);
                S4=case8(S4,N,NM,N1,N2,N3,M,R3,R34,R1,GL3,GLM34,GL1,YLM334,YLM112,Z3,Z3L,Z34,Z34L,Z1,Z1L,F,DF,SDEL,MIN);

                end
            end
        end
    end
end
end

function S4=case1(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,F,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);
    
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
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                for I4=1:2
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                        2*F(I1)*N*valeq*...
                        sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                        *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end 
        end
    end
    
    S4=real(S4);
end

function S4=case2(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,Z1,Z1L,F,DF,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);
    ZE1=Z1(N1);
    ZE1L=Z1L(N1);
    
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
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      2*valeq*(F(I1)*F(I2)*valne1+...
                         F(I1)*DF(I1)*DF(I2)*(1-F(I1))*valne2)*...
                        sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                        *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end
        end
    end
end

function S4=case3(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,Z12,Z12L,F,DF,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);
    ZE12=Z12(N1);
    ZE12L=Z12L(N1);

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
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      2*valeq*(F(I2)*F(I3)*valne1+...
                         F(I2)*DF(I2)*DF(I3)*(1-F(I2))*valne2)*...
                        sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                        *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end
        end
    end
end

function S4=case4(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,Z3,Z3L,F,DF,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);
    ZE3=Z3(N1);
    ZE3L=Z3L(N1);
    
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
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      2*valeq*(F(I3)*F(I4)*valne1+...
                         F(I3)*DF(I3)*DF(I4)*(1-F(I3))*valne2)*...
                        sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                        *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end
        end
    end
end

function S4=case5(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,Z1,Z1L,Z12,Z12L,F,DF,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);            
    ZE1=[Z1(N1),Z1L(N1)];
    ZE12=[Z12(N2),Z12L(N2)];

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
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      2*valeq*(F(I1)*F(I2)*F(I3)*valne(1,1)+...
                          F(I1)*F(I2)*DF(I2)*DF(I3)*(1-F(I2))*valne(1,2)+...
                          F(I1)*DF(I1)*DF(I2)*(1-F(I1))*F(I3)*valne(2,1)+...
                          F(I1)*DF(I1)*DF(I2)*(1-F(I1))*DF(I2)*DF(I3)*(1-F(I2))*valne(2,2))*...
                         sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                         *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end
        end
    end
end

function S4=case6(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,Z12,Z12L,Z3,Z3L,F,DF,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);            
    ZE12=[Z12(N2),Z12L(N2)];
    ZE3=[Z3(N2),Z3L(N2)];

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
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      2*valeq*(F(I2)*F(I3)*F(I4)*valne(1,1)+...
                          F(I2)*F(I3)*DF(I3)*DF(I4)*(1-F(I3))*valne(1,2)+...
                          F(I2)*DF(I2)*DF(I3)*(1-F(I2))*F(I4)*valne(2,1)+...
                          F(I2)*DF(I2)*DF(I3)*(1-F(I2))*DF(I3)*DF(I4)*(1-F(I3))*valne(2,2))*...
                         sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                         *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end
        end
    end
end

function S4=case7(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,Z1,Z1L,Z3,Z3L,F,DF,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);            
    ZE1=[Z1(N2),Z1L(N2)];
    ZE3=[Z3(N2),Z3L(N2)];

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
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      2*valeq*(F(I1)*F(I2)*F(I4)*valne(1,1)+...
                          F(I1)*F(I2)*DF(I2)*DF(I4)*(1-F(I2))*valne(1,2)+...
                          F(I1)*DF(I1)*DF(I2)*(1-F(I1))*F(I4)*valne(2,1)+...
                          F(I1)*DF(I1)*DF(I2)*(1-F(I1))*DF(I2)*DF(I4)*(1-F(I2))*valne(2,2))*...
                         sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                         *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end
        end
    end
end

function S4=case8(S4,N,NM,N1,N2,N3,M,R1,R12,R3,GL1,GLM12,GL3,YLM112,YLM443,Z1,Z1L,Z12,Z12L,Z3,Z3L,F,DF,SDEL,MIN)
    E1=R1(N1);
    E12=R12(N2);
    E3=R3(N3);
    ZE1=[Z1(N2),Z1L(N2)];
    ZE12=[Z12(N2),Z12L(N2)];
    ZE3=[Z3(N3),Z3L(N3)];

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
                token=1;
                % subcase 2
                elseif (abs(ZT1-1)<MIN && abs(ZT12-1)<MIN && abs(ZT3-1)>=MIN)
                    valne(I1,I2,I3)=(-1/6).*((-1)+ZT3).^(-4).*ZT3.*(6+(-11).*N+6.*N.^2+(-1).*N.^3+18.* ...
                           N.*ZT3+(-15).*N.^2.*ZT3+3.*N.^3.*ZT3+(-9).*N.*ZT3.^2+12.*N.^2.* ...
                           ZT3.^2+(-3).*N.^3.*ZT3.^2+2.*N.*ZT3.^3+(-3).*N.^2.*ZT3.^3+N.^3.* ...
                           ZT3.^3+(-6).*ZT3.^N);
                token=2;
                % subcase 3
                elseif (abs(ZT1-1)>=MIN && abs(ZT12-1)<MIN && abs(ZT3-1)<MIN)
                    valne(I1,I2,I3)=(-1/6).*((-1)+ZT1).^(-4).*ZT1.*(12.*ZT1+(-4).*N.*ZT1+(-3).*N.^2.* ...
                          ZT1+N.^3.*ZT1+(-6).*ZT1.^2+2.*N.*ZT1.^2+6.*N.^2.*ZT1.^2+(-2).* ...
                         N.^3.*ZT1.^2+2.*N.*ZT1.^3+(-3).*N.^2.*ZT1.^3+N.^3.*ZT1.^3+(-6).* ...
                          ZT1.^N);
                token=3;
                % subcase 4
                elseif (abs(ZT1-1)<MIN && abs(ZT12-1)>=MIN && abs(ZT3-1)<MIN)
                    valne(I1,I2,I3)=(-1/6).*((-1)+ZT12).^(-4).*ZT12.*(6+(-11).*N+6.*N.^2+(-1).*N.^3+ ...
                   18.*N.*ZT12+(-15).*N.^2.*ZT12+3.*N.^3.*ZT12+(-9).*N.*ZT12.^2+12.* ...
                   N.^2.*ZT12.^2+(-3).*N.^3.*ZT12.^2+2.*N.*ZT12.^3+(-3).*N.^2.* ...
                   ZT12.^3+N.^3.*ZT12.^3+(-6).*ZT12.^N);
                token=4;
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
                token=5;
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
                token=6;
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
                token=7;
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
                token=8;
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
                    S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+SDEL(I1,I2,I3,I4)*...
                      2*valeq*(F(I1)*F(I2)*F(I3)*F(I4)*valne(1,1,1)+...
                                   F(I1)*F(I2)*F(I3)*DF(I3)*DF(I4)*(1-F(I3))*valne(1,1,2)+...
                                   F(I1)*F(I2)*DF(I2)*DF(I3)*(1-F(I2))*F(I4)*valne(1,2,1)+...
                                   F(I1)*F(I2)*DF(I2)*DF(I3)*(1-F(I2))*DF(I3)*DF(I4)*(1-F(I3))*valne(1,2,2)+...
                                   F(I1)*DF(I1)*DF(I2)*(1-F(I1))*F(I3)*F(I4)*valne(2,1,1)+...
                                   F(I1)*DF(I1)*DF(I2)*(1-F(I1))*F(I3)*DF(I3)*DF(I4)*(1-F(I3))*valne(2,1,2)+...
                                   F(I1)*DF(I1)*DF(I2)*(1-F(I1))*DF(I2)*DF(I3)*(1-F(I2))*F(I4)*valne(2,2,1)+...
                                   F(I1)*DF(I1)*DF(I2)*(1-F(I1))*DF(I2)*DF(I3)*(1-F(I2))*DF(I3)*DF(I4)*(1-F(I3))*valne(2,2,2))*...
                         sum(sum((conj(YLM112(:,M+1))*transpose(YLM443(:,M+1)))....
                         *(GL1(:,N1)*transpose(GL3(:,N3))).*GLM12(:,:,M+1,N2)));
                end
            end
        end
    end
end
            
function xf=finite(x)
    MIN=-1e10;
    MAX=1e10;
    xf=x;
    xf(xf>MAX | xf<MIN)=0;
    xf(isnan(xf))=0;
end

function Eig=MatRoots(k,d,ORD)

% find roots of denominator (eigenvalues) by solving eigenvalue problem

if k>8000

    % use large k asmyptotic expansion for large k
    m=(d-3)/2;
    NumPoles=ORD;
    lambda=k;
    alpha=1/sqrt(8*lambda);
    I=complex(0,1);
    r=0;
   for np=1:NumPoles
        if r==NumPoles break; end
        r=r+1;
        n=2*(r-1)+m+1;
        Eig(r)=-(m*(m+1)*0-lambda*I+Epsilon(r-1,d,alpha));
        if imag(Eig(r))~=0
            r=r+1;
            Eig(r)=conj(Eig(r-1));
       end
    end
else

    % use matrix method for intermediate and small k regime
    n=4*ORD;
    E=zeros(n,n);
    for m=1:n
        if k<=1
            a=complex(0,-k*sqrt(m*(m+d-3)/(2*m+d-2)/(2*m+d-4)));            
            if m>1 
                b=complex(0,-k*sqrt((m-1)*((m-1)+d-3)/(2*(m-1)+d-2)/(2*(m-1)+d-4)));
            end
            if m==1
                E(m,1:2)=[(m-1)*(m+d-3),a];
            elseif m==n
                E(m,n-1:n)=[b,(m-1)*(m+d-3)];
            else
                E(m,m-1:m+1)=[b,(m-1)*(m+d-3),a];
            end
        else
            a=complex(0,-sqrt(m*(m+d-3)/(2*m+d-2)/(2*m+d-4)));            
            if m>1 
                b=complex(0,-sqrt((m-1)*((m-1)+d-3)/(2*(m-1)+d-2)/(2*(m-1)+d-4)));
            end
            if m==1
                E(m,1:2)=[(m-1)*(m+d-3)/k,a];
            elseif m==n
                E(m,n-1:n)=[b,(m-1)*(m+d-3)/k];
            else
                E(m,m-1:m+1)=[b,(m-1)*(m+d-3)/k,a];
            end
        end
    end
    TempMat=eig(E);
    [~,index]=sort(real(TempMat));
    TempMat=TempMat(index);
    if k<=1
        Eig=-TempMat(1:ORD);
    else
        Eig=-TempMat(1:ORD)*k;
    end
end
end

function Eig=MatRootsLM(k,N)

Eig=zeros(N,N);

for M=0:(N-1)

% use matrix method for intermediate and small k regime
n=4*N;
E=zeros(n,n);
for I=1:n
    L=I+M;
    if k<=1
        a=complex(0,-k*sqrt((L-M)*(L+M))/sqrt(4*L^2-1));
        if I>1 
		b=complex(0,-k*sqrt((L-1-M)*(L-1+M))/sqrt(4*(L-1)^2-1)); 
	end
        if I==1
            E(I,1:2)=[L*(L-1),a];
        elseif I==n
            E(I,n-1:n)=[b,L*(L-1)];
        else
            E(I,I-1:I+1)=[b,L*(L-1),a];
        end
    else
        a=complex(0,-sqrt((L-M)*(L+M))/sqrt(4*L^2-1));
        if I>1 
		b=complex(0,-sqrt((L-1-M)*(L-1+M))/sqrt(4*(L-1)^2-1)); 
	end
        if I==1
            E(I,1:2)=[L*(L-1)/k,a];
        elseif I==n
            E(I,n-1:n)=[b,L*(L-1)/k];
        else
            E(I,I-1:I+1)=[b,L*(L-1)/k,a];
        end
    end
end
TempMat=eig(E);
[junk,index]=sort(real(TempMat));
TempMat=TempMat(index);
if k<=1
    Eig(:,M+1)=-TempMat(1:N);
else
    Eig(:,M+1)=-TempMat(1:N)*k;
end

end
end

function GLK=gl(K,EigK,ORDEig,ORDL,NumLayer)

GLK=zeros(ORDL,ORDEig);
AL=zeros(NumLayer,1);

for iEig=1:ORDEig
    
    WP=zeros(NumLayer,1);
    dJp=zeros(NumLayer,1);
    WPPROD=zeros(NumLayer,1);
    
    n=NumLayer-1;
    AL(NumLayer)=n/sqrt(4*n^2-1);
    WP(NumLayer)=1/(EigK(iEig)+n*(n+1));
    dJp(NumLayer)=1;
    WPPROD(NumLayer)=1i*K*AL(NumLayer)*WP(NumLayer);

    for n=NumLayer-2:-1:0    
        PL=EigK(iEig)+n*(n+1);        
        AL(n+1)=n/sqrt(4*n^2-1);
        
        WP(n+1)=1/(PL+WP(n+2)*(AL(n+2)*K)^2);
        dJp(n+1)=1-dJp(n+2)*(AL(n+2)*K*WP(n+2))^2;
        WPPROD(n+1)=1i*K*AL(n+1)*WP(n+1);
    end

    GLK(1,iEig)=1/dJp(1);     % L=0 term

    for L=1:(ORDL-1)
        GLK(L+1,iEig)=prod(WPPROD(2:(L+1)))/dJp(1);   % L >= 1
    end
end
end

function [val]=legendrep(P,ORDL)

val=zeros(length(P),ORDL);

val(:,1)=ones(length(P),1);
if ORDL>=2
   val(:,2)=P;
end
   
for N=3:ORDL
    L=N-2;
    val(:,L+1+1)=((2*L+1)*P.*val(:,L+1)-L*val(:,L-1+1))/(L+1);
end

val=transpose(val);
end

function GLMK=glm(K,EigK,ORDEig,ORDL,NumLayer)

GLMK=zeros(ORDL,ORDL,ORDL,ORDEig);
MIN=1e-10;

if abs(K)<MIN
    for M=0:(ORDL-1)
    for L1=M:(ORDL-1)
        GLMK(L1+1,L1+1,M+1,:)=1;
    end
    end

else

for M=0:(ORDL-1)
    AL=zeros(NumLayer,1);
    for iEig=1:ORDEig
        
        WP=zeros(NumLayer,1);
        WM=zeros(NumLayer,1);
        dJp=zeros(NumLayer,1);
        dJm=zeros(NumLayer,1);
        WPPROD=zeros(NumLayer,1);
        
        LP=NumLayer-1;
        AL(NumLayer)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
        WP(NumLayer)=1/(EigK(iEig)+LP*(LP+1));
        
        LM=M;
        if EigK(iEig)+LM*(LM+1)==0
            WM(M+1)=0;
        else
            WM(M+1)=1/(EigK(iEig)+LM*(LM+1));
        end
        AL(M+1)=sqrt((LM-M)*(LM+M))/sqrt(4*LM^2-1);
        
        dJp(NumLayer)=1;
        dJm(M+1)=1;
        WPPROD(NumLayer)=1i*K*AL(NumLayer)*WP(NumLayer);
        
        for n=(NumLayer-2):-1:M
            IP=n+1;
            LP=IP-1;
            IM=NumLayer-n+M;
            LM=IM-1;
            
            PL=EigK(iEig)+LP*(LP+1);
            AL(IP)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
            WP(IP)=1/(PL+WP(IP+1)*(AL(IP+1)*K)^2);
            dJp(IP)=1-dJp(IP+1)*(AL(IP+1)*K*WP(IP+1))^2;
            
            WPPROD(IP)=1i*K*AL(IP)*WP(IP);
            
            PL=EigK(iEig)+LM*(LM+1);
            AL(IM)=sqrt((LM-M)*(LM+M))/sqrt(4*LM^2-1);
            if WM(IM-1)==0
                WM(IM)=0;
            else
                WM(IM)=1/(PL+WM(IM-1)*(AL(IM)*K)^2);
            end
            dJm(IM)=1-dJm(IM-1)*(AL(IM)*K*WM(IM-1))^2;
            
        end
        
        for L1=M:(ORDL-1)
            for L2=L1:(ORDL-1)
                L=sort([L1 L2]);
                IM=L(1)+1;
                IP=L(2)+1;
                if L1==M
                    dWL0M=1/dJp(IM);
                else
                    dWL0M=1/(-(AL(IM)*K*WM(IM-1))^2*dJm(IM-1)+dJp(IM));
                end

                if IM==IP
                    GLMK(L1+1,L2+1,M+1,iEig)=dWL0M;
                else
                    GLMK(L1+1,L2+1,M+1,iEig)=dWL0M*prod(WPPROD((IM+1):IP));
                end
                GLMK(L2+1,L1+1,M+1,iEig)=GLMK(L1+1,L2+1,M+1,iEig);
            end
        end
        
    end
end
end
end

function [RHO,PHI]=rotvec(E1,E2)
% E1 and E2 are unit vectors
B1=acos(E1(3));
A1=atan2(E1(2),E1(1)); 

ROTZ=[cos(-A1) -sin(-A1) 0;sin(-A1) cos(-A1) 0;0 0 1];
ROTY=[cos(-B1) 0 sin(-B1);0 1 0;-sin(-B1) 0 cos(-B1)];

E2=ROTY*ROTZ*E2;

RHO=E2(3);
PHI=atan2(E2(2),E2(1));
end

