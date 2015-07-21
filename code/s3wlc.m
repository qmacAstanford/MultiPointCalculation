function S3=s3wlc(N,NM,LAM,FA,Q1,Q2,Q3,ORDEig,ORDL,NumLayer)

S3=zeros(2,2,2);
MIN=5e-4;
digits(5);

if sum(power(Q1+Q2+Q3,2)) <= MIN
    
    % Evaluate the quantities for s3 calculation
    FB=1-FA;
    F=[FA,FB];
    DF=[1,-1];
    
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));

    EQ1=Q1/Q1MAG;
    EQ2=Q2/Q2MAG;
    EQ3=Q3/Q3MAG;
    
    RHO12=sum(EQ1.*EQ2);
    RHO23=sum(EQ2.*EQ3);
    RHO13=sum(EQ1.*EQ3);
    
    PL12=legendrep(-RHO12,ORDL);
    PL23=legendrep(-RHO23,ORDL);
    PL13=legendrep(-RHO13,ORDL);
   
    R1=MatRoots(Q1MAG,3,ORDEig);
    R2=MatRoots(Q2MAG,3,ORDEig);
    R3=MatRoots(Q3MAG,3,ORDEig);
    
    GL1=gl(Q1MAG,R1,ORDEig,ORDL,NumLayer);
    GL2=gl(Q2MAG,R2,ORDEig,ORDL,NumLayer);
    GL3=gl(Q3MAG,R3,ORDEig,ORDL,NumLayer);
    
    Z1=exp(R1*NM);
    Z2=exp(R2*NM);
    Z3=exp(R3*NM);
    
    Z1L=exp(R1*NM)*LAM;
    Z2L=exp(R2*NM)*LAM;
    Z3L=exp(R3*NM)*LAM;
    
    for N1=1:ORDEig
        for N2=1:ORDEig
            
    % Case 1: J1=J2=J3
        S3=case1(S3,N,NM,N1,N2,R1,R2,Z1,Z2,PL12,GL1,GL2,F,MIN);
        S3=case1(S3,N,NM,N1,N2,R1,R3,Z1,Z3,PL13,GL1,GL3,F,MIN);
        S3=case1(S3,N,NM,N1,N2,R2,R3,Z2,Z3,PL23,GL2,GL3,F,MIN);

    % Case 2: J1<J2=J3
        % J1<J2=J3
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(2,1,1)=1;
        SDEL(1,2,2)=1;
        S3=case2(S3,N,NM,N1,N2,R1,R2,Z1,Z2,Z1L,PL12,GL1,GL2,F,DF,SDEL,1,2,MIN);
        S3=case2(S3,N,NM,N1,N2,R1,R3,Z1,Z3,Z1L,PL13,GL1,GL3,F,DF,SDEL,1,2,MIN);

        % J2<J1=J3
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(2,1,2)=1;
        SDEL(1,2,1)=1;
        S3=case2(S3,N,NM,N1,N2,R2,R3,Z2,Z3,Z2L,PL23,GL2,GL3,F,DF,SDEL,1,2,MIN);
        S3=case2(S3,N,NM,N1,N2,R2,R1,Z2,Z1,Z2L,PL12,GL2,GL1,F,DF,SDEL,1,2,MIN);

        % J3<J1=J2
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,1,2)=1;
        SDEL(2,2,1)=1;
        S3=case2(S3,N,NM,N1,N2,R3,R1,Z3,Z1,Z3L,PL13,GL3,GL1,F,DF,SDEL,3,1,MIN);
        S3=case2(S3,N,NM,N1,N2,R3,R2,Z3,Z2,Z3L,PL23,GL3,GL2,F,DF,SDEL,3,1,MIN);

    % Case 3: J1=J2<J3
        % J1=J2<J3
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,1,2)=1;
        SDEL(2,2,1)=1;
        S3=case3(S3,N,NM,N1,N2,R1,R3,Z1,Z3,Z3L,PL13,GL1,GL3,F,DF,SDEL,2,3,MIN);
        S3=case3(S3,N,NM,N1,N2,R2,R3,Z2,Z3,Z3L,PL23,GL2,GL3,F,DF,SDEL,2,3,MIN);

        % J1=J3<J2
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,2,1)=1;
        SDEL(2,1,2)=1;
        S3=case3(S3,N,NM,N1,N2,R1,R2,Z1,Z2,Z2L,PL12,GL1,GL2,F,DF,SDEL,3,2,MIN);
        S3=case3(S3,N,NM,N1,N2,R3,R2,Z3,Z2,Z2L,PL23,GL3,GL2,F,DF,SDEL,3,2,MIN);

        % J2=J3<J1
        SDEL=zeros(2,2,2);
        SDEL(1,1,1)=1;
        SDEL(2,2,2)=1;
        SDEL(1,2,2)=1;
        SDEL(2,1,1)=1;
        S3=case3(S3,N,NM,N1,N2,R2,R1,Z2,Z1,Z1L,PL12,GL1,GL2,F,DF,SDEL,3,1,MIN);
        S3=case3(S3,N,NM,N1,N2,R3,R1,Z3,Z1,Z1L,PL23,GL2,GL3,F,DF,SDEL,3,1,MIN);

    % Case 4: J1<J2<J3
        SDEL=ones(2,2,2);
        % J1<J2<J3
        S3=case4(S3,N,NM,N1,N2,R1,R3,Z1,Z3,Z1L,Z3L,PL13,GL1,GL3,F,DF,SDEL,1,2,3,MIN);
        
        %J3<J2<J1
        S3=case4(S3,N,NM,N1,N2,R3,R1,Z3,Z1,Z3L,Z1L,PL13,GL3,GL1,F,DF,SDEL,3,2,1,MIN);
        
        %J1<J3<J2
        S3=case4(S3,N,NM,N1,N2,R1,R2,Z1,Z2,Z1L,Z2L,PL12,GL1,GL2,F,DF,SDEL,1,3,2,MIN);
        
        %J2<J3<J1
        S3=case4(S3,N,NM,N1,N2,R2,R1,Z2,Z1,Z2L,Z1L,PL12,GL2,GL1,F,DF,SDEL,2,3,1,MIN);
        
        %J2<J1<J3
        S3=case4(S3,N,NM,N1,N2,R2,R3,Z2,Z3,Z2L,Z3L,PL23,GL2,GL3,F,DF,SDEL,2,1,3,MIN);
        
        %J3<J1<J2
        S3=case4(S3,N,NM,N1,N2,R3,R2,Z3,Z2,Z3L,Z2L,PL23,GL3,GL2,F,DF,SDEL,3,1,2,MIN);
        end
    end
    S3=real(S3);
end
end

function S3=case1(S3,N,NM,N1,N2,E1,E2,ZT1,ZT2,PL12,GL1,GL2,F,MIN)
% Case 1: J1=J2=J3
    R1=E1(N1);
    R2=E2(N2);
    Z1=ZT1(N1);
    Z2=ZT2(N2);
    
    % on same monomer
    if abs(R1-R2)<MIN
        valeq=R1.^(-3).*(2+NM.*R1+((-2)+NM.*R1).*Z1);
    else
        valeq=R1.^(-2).*(R1+(-1).*R2).^(-1).*R2.^(-2).*((-1).*NM.*R1.*R2.^2+ ...
            R2.^2.*((-1)+Z1)+R1.^2.*(1+NM.*R2+(-1).*Z2));
    end
    
    S3(1,1,1)=S3(1,1,1)+2*F(1)*N*valeq*sum(PL12.*GL1(:,N1).*GL2(:,N2));
    S3(2,2,2)=S3(2,2,2)+2*F(2)*N*valeq*sum(PL12.*GL1(:,N1).*GL2(:,N2));
end

function S3=case2(S3,N,NM,N1,N2,E1,E2,ZT1,ZT2,ZT1L,PL12,GL1,GL2,F,DF,SDEL,A1,A2,MIN)
% Case 2: J1<J2=J3
    R1=E1(N1);
    R2=E2(N2);
    Z1=ZT1(N1);
    Z2=ZT2(N2);
    Z1L=ZT1L(N1);
    
    % on same monomer
    if abs(R1-R2)<MIN
        valeq=R1.^(-3).*((-1)+Z1).*Z1.^(-1).*(1+((-1)+NM.*R1).*Z1);
    else
        valeq=R1.^(-2).*(R1+(-1).*R2).^(-1).*R2.^(-1).*((-1)+Z1).*Z1.^(-1).*(R1+ ...
            R2.*((-1)+Z1)+(-1).*R1.*Z2);
    end
    
    % on different monomers
    if abs(Z1-1)<MIN
        valne1=(1/2).*((-1).*N+N.^2);
    else
        valne1=(-1).*((-1)+Z1).^(-2).*Z1.*(1+N.*((-1)+Z1)+(-1).*Z1.^N);
    end
    
    if abs(Z1L-1)<MIN
        valne2=(1/2).*((-1).*N+N.^2);
    else
        valne2=(-1).*((-1)+Z1L).^(-2).*Z1L.*(1+N.*((-1)+Z1L)+(-1).*Z1L.^N);
    end
    valeq=finite(valeq);
    valne1=finite(valne1);
    valne2=finite(valne2);
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                          valeq*(FV(A1)*FV(A2)*valne1+...
                             FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*valne2)*...
                          sum(PL12.*GL1(:,N1).*GL2(:,N2));
            end
        end
    end
end

function S3=case3(S3,N,NM,N1,N2,E1,E2,ZT1,ZT2,ZT2L,PL12,GL1,GL2,F,DF,SDEL,A2,A3,MIN)
% Case 3: J1=J2<J3
    R1=E1(N1);
    R2=E2(N2);
    Z1=ZT1(N1);
    Z2=ZT2(N2);
    Z2L=ZT2L(N2);
    
    % on same monomer
    if abs(R1-R2)<MIN
        valeq=R1.^(-3).*((-1)+Z1).*Z1.^(-1).*(1+((-1)+NM.*R1).*Z1);
    else
        valeq=R1.^(-1).*R2.^(-2).*((-1).*R1+R2).^(-1).*(R2+(-1).*R2.*Z1+R1.*(( ...
            -1)+Z2)).*((-1)+Z2).*Z2.^(-1);
    end
    
    % on different monomers
    if abs(Z2-1)<MIN
        valne1=(1/2).*((-1).*N+N.^2);
    else
        valne1=(-1).*((-1)+Z2).^(-2).*Z2.*(1+N.*((-1)+Z2)+(-1).*Z2.^N);
    end
    if abs(Z2L-1)<MIN
        valne2=(1/2).*((-1).*N+N.^2);
    else
        valne2=(-1).*((-1)+Z2L).^(-2).*Z2L.*(1+N.*((-1)+Z2L)+(-1).*Z2L.^N);
    end
    valeq=finite(valeq);
    valne1=finite(valne1);
    valne2=finite(valne2);
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                          valeq*(FV(A2)*FV(A3)*valne1+...
                             FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne2)*...
                          sum(PL12.*GL1(:,N1).*GL2(:,N2));
            end
        end
    end
end

function S3=case4(S3,N,NM,N1,N2,E1,E2,ZT1,ZT2,ZT1L,ZT2L,PL12,GL1,GL2,F,DF,SDEL,A1,A2,A3,MIN)
% Case 4: J1<J2<J3
    R1=E1(N1);
    R2=E2(N2);
    Z1=ZT1(N1);
    Z2=ZT2(N2);
    Z1L=ZT1L(N1);
    Z2L=ZT2L(N2);
    
    ZE1=[Z1,Z1L];
	ZE2=[Z2,Z2L];

    % on same monomer
    if abs(R1-R2)<MIN
        valeq=2.*NM.*R1.^(-2).*((-1)+(1/2).*(Z1.^(-1)+Z1));
    else
        valeq=R1.^(-1).*(R1+(-1).*R2).^(-1).*R2.^(-1).*((-1)+Z1).*Z1.^(-1).*(Z1+ ...
          (-1).*Z2).*((-1)+Z2).*Z2.^(-1);
    end
    valeq(valeq>1e5)=0;
    
    % on different monomers
%     MIN=abs(exp(NM)*exp(R1)*(1-exp(MIN)));
    valne=zeros(2,2);
    %valne(1,1) = no lambda term
    %valne(1,2) = lambda^(J3-J2) term
    for I1=1:2
        for I2=1:2
            ZT1=ZE1(I1);
            ZT2=ZE2(I2);
            if (abs(ZT1-1)<MIN && abs(ZT2-1)>=MIN)
                valne(I1,I2)=(-1/2).*((-1)+ZT2).^(-3).*ZT2.*(2+N.^2.*((-1)+ZT2).^2+(-2).* ...
                ZT2.^N+(-1).*N.*(3+(-4).*ZT2+ZT2.^2));
            elseif (abs(ZT1-1)>=MIN && abs(ZT2-1)<MIN)
                valne(I1,I2)=(1/2).*(2+(-3).*N+N.^2).*((-1)+ZT1).^(-2).*ZT1+...
                    (-1/2).*((-1)+ZT1).^(-3).*ZT1.*(2.*ZT1+N.*ZT1+(-1).*N.^2.*ZT1+(-1) ...
                    .*N.*ZT1.^2+N.^2.*ZT1.^2+(-2).*ZT1.^N);
            elseif (abs(ZT1-1)<MIN && abs(ZT2-1)<MIN)
                valne(I2,I2)=(1/6).*N.*(2+(-3).*N+N.^2);
            elseif abs(ZT1-ZT2)<MIN
                valne(I1,I2)=((-1)+ZT1).^(-3).*ZT1.*(2.*ZT1+(-1).*N.*ZT1+N.*ZT1.^2+(-1).*N.* ...
                ZT1.^N+(-2).*ZT1.^(1+N)+N.*ZT1.^(1+N));
            else
                valne(I1,I2)=((-1)+ZT1).^(-2).*ZT1.*(ZT1+(-1).*ZT2).^(-1).*((-1)+ZT2).^(-2).* ...
                 ZT2.*(ZT1.^N.*((-1)+ZT2).^2+(-1).*((-2)+N).*ZT2+((-1)+N).*ZT2.^2+( ...
                  -1).*ZT2.^N+(-1).*ZT1.^2.*((-1)+N+(-1).*N.*ZT2+ZT2.^N)+ZT1.*((-2)+ ...
                  N+(-1).*N.*ZT2.^2+2.*ZT2.^N));
            end
        end
    end
    valeq=finite(valeq);
    valne=finite(valne);
    
    for I1=1:2
        for I2=1:2
            for I3=1:2
                FV=[F(I1),F(I2),F(I3)];
                DFV=[DF(I1),DF(I2),DF(I3)];
                S3(I1,I2,I3)=S3(I1,I2,I3)+SDEL(I1,I2,I3)*...
                   valeq*(FV(A1)*FV(A2)*FV(A3)*valne(1,1)+...
                      FV(A1)*FV(A2)*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(1,2)+...
                      FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*FV(A3)*valne(2,1)+...
                      FV(A1)*DFV(A1)*DFV(A2)*(1-FV(A1))*DFV(A2)*DFV(A3)*(1-FV(A2))*valne(2,2))*...
                   sum(PL12.*GL1(:,N1).*GL2(:,N2));
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

Eig=zeros(ORD,1);

if k>8000

    % use large k asmyptotic expansion for large k
    
    for I=1:floor(ORD/2)
        l=I-1;
        alpha=1/sqrt(8*k);
        Eig(2*l+1)=1i*k-Epsilon(l,d,alpha);
        Eig(2*l+2)=conj(Eig(2*l+1));
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

function value=Epsilon(l,d,alpha)

% eigenvalues using large k asymptotic expansion
% generates epsilon^{s}_r (see the paper)

I=complex(0,1);
beta=-sqrt(2)/4*(1+I);
m=(d-3)/2;
n=2*l+m+1;

epsilon_0=(-1/2/beta)^(-1)*(n/2);
epsilon_1=(-1/2/beta)^( 0)*(-1/8*(n^2+3-3*m^2)-m*(m+1));
epsilon_2=(-1/2/beta)^( 1)*(-1/2^5*n*(n^2+3-9*m^2));
epsilon_3=(-1/2/beta)^( 2)*(-1/2^8*(5*n^4+34*n^2+9)-(102*n^2+42)*m^2+33*m^4);
epsilon_4=(-1/2/beta)^( 3)*(-1/2^11*n*(33*n^4+410*n^2+405)-(1230*n^2+1722)*m^2+813*m^4);
epsilon_5=(-1/2/beta)^( 4)*(-1/2^12*9*(7*n^6+140*n^4+327*n^2+54-(420*n^4+1350*n^2+286)*m^2+(495*n^2+314)*m^4-82*m^6));

value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha^2+...
      epsilon_4*alpha^3+epsilon_5*alpha^4;
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
