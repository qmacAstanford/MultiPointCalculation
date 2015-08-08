function S4=s4wlcV2(N,NM,LAM,FA,Q1,Q2,Q3,Q4,ORDEig,ORDL,NumLayer)
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

% exclude zero vector
MIN=1e-6;
if sum(power(Q1+Q2+Q3+Q4,2)) > MIN
    return
end

% Evaluate the quantities for s4 calculation
FB=1-FA;
F=[FA,FB];


S4=zeros(2,2,2,2);
%S4split=zeros(8,2,2,2,2);
orders = perms(1:4);
for orderNum=1:24
    order=orders(orderNum,:);
    Q = [Q1,Q2,Q3,Q4];
    Qnew=[Q(:,order(1)),Q(:,order(2)),Q(:,order(3)),Q(:,order(4))];
    % Qnew is the reordered Q

    % Now take advantage of translatioanl invariance
    q1=-Qnew(:,1);
    q2=-(Qnew(:,1)+Qnew(:,2));
    q3=Qnew(:,4);

    % Find unit vectors
    Q1_n=q1/norm(q1,2);
    Q2_n=q2/norm(q2,2);
    Q3_n=q3/norm(q3,2);

    % Find Euler angles
    aligned=10^-13; % angle in radians betwene to vectors before I assume they are the same
    if (Q2_n-Q1_n) < aligned
        cosB1=0;
        alpha1=0;
        cosB2=dot(Q2_n,Q3_n);
    elseif (Q2_n-Q1_n) < aligned
        cosB1=dot(Q2_n,Q1_n);
        alpha1=0;
        cosB2=0;
    else
        cosB1=dot(Q2_n,Q1_n);
        cosB2=dot(Q2_n,Q3_n);
        v1=cross(Q2_n,Q1_n)/norm(cross(Q2_n,Q1_n));
        v2=cross(Q2_n,Q3_n)/norm(cross(Q2_n,Q3_n));
        alpha1=acos(dot(v1,v2));
    end 
    % Calculate Wigner D matrices
    firstD=WignerD_lm0(L_max,alpha1,acos(cosB1));
    secondD=WingerD_lm0(L_max,0,acos(cosB2));
    
    
    % Now calculate the eigenvalues 
    R1=MatRoots(norm(q1),3,ORDEig);  % 1 x ORDEig  
    R12=MatRootsLM(norm(q2),ORDEig); % ORDEig x ORDEig
    % dimensions (L+1,M+1)
    R4=MatRoots(norm(q4),3,ORDEig);  % 1 x OREig  
%     R1=-dot(q1,q1)/(2*d);
%     R12=-dot(q2,q2)/(2*d);
%     R4=-dot(q3,q3)/(2*d);
    GL1=gl(norm(q1),R1,ORDEig,ORDL,NumLayer);   % ORDL x ORDEig
    GLM12=glm(norm(q2),R12,ORDEig,ORDL,NumLayer); % ORDL x ORDL x ORDL x ORDEig
    % dimensions (L+1,L+1,M+1,N2)
    GL4=gl(norm(q3),R4,ORDEig,ORDL,NumLayer);   % ORDL x ORDEig
    
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

