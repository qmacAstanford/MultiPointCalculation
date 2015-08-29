function S4=s4wlcV2(N,NM,LAM,FA,Q1,Q2,Q3,Q4,ORDEig,ORDL,NumLayer)
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
    if norm(q2) < aligned
        % if q2=0 we may as well point it perpendicular to the other two
        cosB1=0;
        cosB2=0;
        alpha1=dot(Q1_n,Q3_n);
    elseif abs(Q2_n-Q1_n) < aligned
        cosB1=1;
        alpha1=0;
        cosB2=dot(Q2_n,Q3_n);
    elseif abs(Q2_n-Q3_n) < aligned
        cosB1=dot(Q2_n,Q1_n);
        alpha1=0;
        cosB2=1;
    else
        cosB1=dot(Q2_n,Q1_n);
        cosB2=dot(Q2_n,Q3_n);
        v1=cross(Q2_n,Q1_n)/norm(cross(Q2_n,Q1_n));
        v2=cross(Q2_n,Q3_n)/norm(cross(Q2_n,Q3_n));
        alpha1=acos(dot(v1,v2));
    end 
    % Calculate Wigner D matrices
    firstD=WignerD_lm0(ORDL,alpha1,cosB1);
    secondD=WignerD_lm0(ORDL,0,cosB2);
    % index: (lam+1,M+1)
    
    % Now calculate the eigenvalue
    R1mtrx=Eigenvalues(norm(q1),ORDEig,ORDL); % won't need all of output
    R12mtrx=Eigenvalues(norm(q2),ORDEig,ORDL);
    R4mtrx=Eigenvalues(norm(q3),ORDEig,ORDL); % won't need all of output
    % index: (l+1,M+1)
    
    GL1=Residues(norm(q1),R1mtrx,ORDEig,ORDL,NumLayer,1);  % won't need all of output
    GLM12=Residues(norm(q2),R12mtrx,ORDEig,ORDL,NumLayer,ORDL);
    GL4=Residues(norm(q3),R4mtrx,ORDEig,ORDL,NumLayer,1);  % won't need all of output
    % index: (lam1+1, lam2+1, mu+1, l+1)
    
    
    for M=0:(ORDL-1)
        if mod(ORDEig-M,2)==0
            Lmax=ORDEig-1;
        else
            Lmax=ORDEig-2;
        end
        if Lmax < M
            error('make ORDEig larger')
        end
        for L1=M:Lmax
            for L2=M:Lmax
                for L3=M:Lmax
                    
                    R1=  R1mtrx(L1+1,M+1);
                    R12=R12mtrx(L2+1,M+1);
                    R4=  R4mtrx(L3+1,M+1);
                    
                    Z1= exp(R1*NM);
                    Z12=exp(R12*NM);
                    Z4= exp(R4*NM);
                    Z1L=  Z1*LAM;
                    Z12L=Z12*LAM;
                    Z4L=  Z4*LAM;
                    
%                     S4term=zeros(2,2,2,2);
%                     S4term=S4term+case1(N,NM,R1,R12,R4,F);
%                     S4term=S4term+case2(N,NM,R1,R12,R4,Z1,Z1L,                F,order);
%                     S4term=S4term+case3(N,NM,R1,R12,R4,       Z12,Z12L,       F,order);
%                     S4term=S4term+case4(N,NM,R1,R12,R4,                Z4,Z4L,F,order);
%                     S4term=S4term+case5(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,       F,order);
%                     S4term=S4term+case6(N,NM,R1,R12,R4,       Z12,Z12L,Z4,Z4L,F,order);
%                     S4term=S4term+case7(N,NM,R1,R12,R4,Z1,Z1L,         Z4,Z4L,F,order);
%                     S4term=S4term+case8(N,NM,R1,R12,R4,Z1,Z1L,Z12,Z12L,Z4,Z4L,F,order);

                    valsne=zeros(2,2,2,8);
                    valseq=zeros(8,1);
                    
                    valseq(1)=case1Int(R1,R12,R4,NM);
                    valsne(:,:,:,1)=ones(2,2,2)*N;
                    
                    valseq(2)=case4Int(R4,R12,R1,NM);
                    valsne(:,:,:,2)=case2sum(N,Z1,Z1L);
                    
                    valseq(3)=case3Int(R1,R12,R4,NM);
                    valsne(:,:,:,3)=case3sum(N,Z12,Z12L);
                    
                    valseq(4)=case4Int(R1,R12,R4,NM); % Case 4: J1==J2==J3 <J4
                    valsne(:,:,:,4)=case4sum(N,Z4,Z4L);
                    
                    valseq(5)=case6Int(R4,R12,R1,NM); % Case 5: J1 <J2 <J3==J4
                    valsne(:,:,:,5)=case5sum(N,Z1,Z1L,Z12,Z12L);
                    
                    valseq(6)=case6Int(R1,R12,R4,NM); % Case 6: J1==J2 <J3 <J4
                    valsne(:,:,:,6)=case6sum(N,Z12,Z12L,Z4,Z4L);
                    
                    valseq(7)=case7Int(R1,R12,R4,NM);% Case 7: J1 <J2==J3 <J4
                    valsne(:,:,:,7)=case7sum(N,Z1,Z1L,Z4,Z4L);
                    
                    valseq(8)=case8Int(R1,R12,R4,NM);
                    valsne(:,:,:,8)=case8sum(N,Z1,Z1L,Z12,Z12L,Z4,Z4L);
                    
                    S4term=BinomialSum(valseq,valsne,order,F);
                    
                    if isnan(S4term)
                        error('S4 term is NaN')
                    end
                    RotTerm=0;
                    if M==0 % If statment needed to include negitive \mu values
                        for lam2=M:(ORDL-1)
                            for lam3=M:(ORDL-1)
                                RotTerm=RotTerm+firstD(lam2+1,M+1)...
                                               *conj(secondD(lam3+1,M+1))...
                                               *GL1(1,lam2+1,M+1,L1+1)...
                                               *GLM12(lam2+1,lam3+1,M+1,L2+1)...
                                               *GL4(lam3+1,1,M+1,L3+1);
                                if isnan(RotTerm)
                                    sprintf('firstD(lam2+1,M+1)=%g+%gi',real(firstD(lam2+1,M+1)),imag(firstD(lam2+1,M+1)))
                                    sprintf('secondD(lam3+1,M+1)=%g+%gi',real(secondD(lam3+1,M+1)),imag(secondD(lam3+1,M+1)))
                                    sprintf('GL1(1,lam2+1,M+1,L1+1)=%g+%gi',real(GL1(1,lam2+1,M+1,L1+1)),imag(GL1(1,lam2+1,M+1,L1+1)))
                                    sprintf('GLM12(lam2+1,lam3+1,M+1,L2+1)=%g+%gi',real(GLM12(lam2+1,lam3+1,M+1,L2+1)),imag(GLM12(lam2+1,lam3+1,M+1,L2+1)))
                                    sprintf('GL4(lam3+1,1,M+1,L3+1)=%g+%gi',real(GL4(lam3+1,1,M+1,L3+1)),imag(GL4(lam3+1,1,M+1,L3+1)))
                                    error('NaN encountered')
                                end
                            end
                        end
                    else
                        for lam2=M:(ORDL-1)
                            for lam3=M:(ORDL-1)
                                RotTerm=RotTerm+2*real(firstD(lam2+1,M+1)...
                                               *conj(secondD(lam3+1,M+1)) )...
                                               *GL1(1,lam2+1,1,L1+1)...
                                               *GLM12(lam2+1,lam3+1,M+1,L2+1)...
                                               *GL4(lam3+1,1,1,L3+1);
                                           
                                if isnan(RotTerm)
                                    sprintf('lam2=%d, lam3=%d, M=%d',lam2,lam3,M)
                                    sprintf('firstD(lam2+1,M+1)=%g+%gi',real(firstD(lam2+1,M+1)),imag(firstD(lam2+1,M+1)))
                                    sprintf('secondD(lam3+1,M+1)=%g+%gi',real(secondD(lam3+1,M+1)),imag(secondD(lam3+1,M+1)))
                                    sprintf('GL1(1,lam2+1,M+1,L1+1)=%g+%gi',real(GL1(1,lam2+1,M+1,L1+1)),imag(GL1(1,lam2+1,M+1,L1+1)))
                                    sprintf('GLM12(lam2+1,lam3+1,M+1,L2+1)=%g+%gi',real(GLM12(lam2+1,lam3+1,M+1,L2+1)),imag(GLM12(lam2+1,lam3+1,M+1,L2+1)))
                                    sprintf('GL4(lam3+1,1,M+1,L3+1)=%g+%gi',real(GL4(lam3+1,1,M+1,L3+1)),imag(GL4(lam3+1,1,M+1,L3+1)))
                                    error('NaN encountered')
                                end
                            end
                        end
                    end
                    
%                     v2= ( squeeze(GL1(1,(M+1):ORDL,1,L1))'.*firstD((M+1):ORDL,M+1) )'; % row mtrx by lam2
%                     v3= ( squeeze(GL4((M+1):ORDL,1,1,L3)).*conj(secondD((M+1):ORDL,M+1)) ); % column mtrx by lam3
%                     mtrx= squeeze(GLM12((M+1):ORDL,(M+1):ORDL,M+1,L2+1); % matrix (lam2,lam3)
%                     RotTerm=v2*mtrx*v3;                   
                    S4=S4+sum(S4term,5)*RotTerm;
%                     sprintf('M=%d,K1=%g,K2=%g,K3=%g, L1=%d, L2=%d, L3=%d, S4term=%g+%gi, RotTerm=%g+%gi',...
%                             M,norm(q1),norm(q2),norm(q3),L1,L2,L3,real(S4term(1,1,1,1)),imag(S4term(1,1,1,1)),real(RotTerm),imag(RotTerm))
                end
            end
        end
    end
    
   
end

end 

function Eig=Eigenvalues(K,ORDEig,ORDL)
% Calculates the eigenvalues for the inverse laplace transform
% Returns a ORDEig x ORDL matrix
% Index of output: (l+1,M+1)
% Choose cutoff manually based on previous compairison
% You may not need all output values, i.e. if you know M=0

cutoff = 800; % This appears to be a pretty good cutoff

if K<cutoff
    Eig=IntermediateKEigenValues(K,ORDEig,ORDL);
else
    Eig=LargeKEigenValues(K,ORDEig,ORDL,3);
end

end

function Eig=IntermediateKEigenValues(k,ORDEig,ORDL)
% essentially does the same thing as MatRootsLM
% output size: ORDEig x ORDL
% output indices: (l+1,m+1)

N=2*(ceil(ORDEig/2));
Eig=zeros(N,ORDL);  % output no-longer square

for M=0:(ORDL-1)
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
    if N-M<1
        Eig(:,M+1)=zeros(N,1)*NaN;
    elseif k<=1
        Eig(:,M+1)=-TempMat(1:N); 
    else
        Eig(:,M+1)=-TempMat(1:N)*k; 
    end
end

% Inforce choice of how to order l values, choice is somewhat arbitrary
% This used to be done outside of MatRoots
for M=0:ORDL-1
    for I=1:2:ORDEig
        Eig(I,M+1)=real(Eig(I,M+1))+1i*abs(imag(Eig(I,M+1)));
        Eig(I+1,M+1)=real(Eig(I+1,M+1))-1i*abs(imag(Eig(I+1,M+1)));    
    end
end

% delete extra row if it exists
if N>ORDEig
    Eig(end,:)=[];
end

% shift to higher L
for M=0:(ORDL-1)
    if ORDEig-M<1
        Eig(:,M+1)=zeros(ORDEig,1)*NaN;
    else  
        Eig(:,M+1)=[zeros(M,1)*NaN;Eig(1:(ORDEig-M),M+1)]; % Q.J.M. added this 8/16/15
    end
end


end

function Eig=LargeKEigenValues(k,ORDEig,ORDL,d)
% Returs an ORDEig x ORDL  where entries are for l+1,mu+1
% where l is the order of the eigenvalue
% and mu is is the azimuthal quantum number (i.e. m)

Eig=zeros(ORDEig,ORDL);
mu=0:(ORDL-1);

for I=1:floor(ORDEig/2)
    l=I-1;
    alpha=1/sqrt(8*k);
    Eig(2*l+1,:)=1i*k-mu.*(mu+d-2)-Epsilon(l,d,alpha,mu); % Q.J.M. changes this line 8/1/15
    Eig(2*l+2,:)=conj(Eig(2*l+1,:));
end

% Q.J.M. added this loop 8/16/15
for M=0:(ORDL-1)
    if ORDEig-M<1
        Eig(:,M+1)=zeros(ORDEig,1)*NaN;
    else  
        Eig(:,M+1)=[zeros(M,1)*NaN;Eig(1:(ORDEig-M),M+1)];
    end
end
    
end

function value=Epsilon(l,d,alpha,mu)
% For use in calculating Large K EigenValues
% eigenvalues using large k asymptotic expansion
% generates epsilon^{s}_r (see the paper)

I=complex(0,1);
beta=-sqrt(2)/4*(1+I);
m=mu+(d-3)/2;  % Q.J.M. changes this line 8/1/15
n=2*l+m+1;  % also know as s

epsilon_0=(-1/2/beta)^(-1)*(n/2);
epsilon_1=(-1/2/beta)^( 0)*(-1/8*(n.^2+3-3*m.^2)-m.*(m+1));
epsilon_2=(-1/2/beta)^( 1)*(-1/2^5*n.*(n.^2+3-9*m.^2));
epsilon_3=(-1/2/beta)^( 2)*(-1/2^8*(5*n.^4+34*n.^2+9)-(102*n.^2+42).*m.^2+33*m.^4);
epsilon_4=(-1/2/beta)^( 3)*(-1/2^11*n.*(33*n.^4+410*n.^2+405)-(1230*n.^2+1722).*m.^2+813*m.^4);
epsilon_5=(-1/2/beta)^( 4)*(-1/2^12*9*(7*n.^6+140*n.^4+327*n.^2+54-(420*n.^4+1350*n.^2+286).*m.^2+(495*n.^2+314).*m.^4-82*m.^6));

value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha^2+...
      epsilon_4*alpha^3+epsilon_5*alpha^4;
end

function out=WignerD_lm0(ORDL,alpha,cosB)
% Wigner D Matrix when second L index is zero
% Returns a ORDL x ORDL matrix with indicies (lam+1,M+1)

out=zeros(ORDL,ORDL)*NaN;

for L=0:(ORDL-1)
    M=(0:L);
    out(L+1,1:(L+1))=legendre(L,cosB)'.*...
                     sqrt(factorial(L-M)./factorial(L+M))...
                     .*exp(1i*M*alpha);
end
end

function out=Residues(K,EigK,ORDEig,ORDL,NumLayer,ORDMU)
% To speed this function up you may not need to return all M values,
% you could separate out a ORDM
% Returns a ORDL x ORDL x ORDMU x ORDEig matrix
% Index (lam1+1, lam2+1, mu+1, l+1)
% I haved included the ORDMU separate from ORDL to save computation time

cutoff=10^(-11); % chosen by looking at graph

largeK=glm(K,EigK,ORDEig,ORDL,NumLayer,ORDMU);

smallK=zeros(ORDL,ORDL,ORDL,ORDEig)*NaN;
out=zeros(ORDL,ORDL,ORDMU,ORDEig)*NaN;

for lam1=0:ORDL-1
    for lam2=0:ORDL-1
        for iEig=1:ORDEig
            for mu=0:min([lam1,lam2,ORDMU-1])

                smallK(lam1+1,lam2+1,mu+1,iEig)=SmallAsympRes(K,iEig,lam1,lam2,mu,3);
                if abs(smallK(lam1+1,lam2+1,mu+1,iEig)) < cutoff
                    out(lam1+1,lam2+1,mu+1,iEig)=smallK(lam1+1,lam2+1,mu+1,iEig);
                else
                    out(lam1+1,lam2+1,mu+1,iEig)=largeK(lam1+1,lam2+1,mu+1,iEig);
                end
            end
        end
    end
end


end

function GLK=SmallAsympRes(K,iEig,lam1,lam2,mu,d)
% calculate the residual using small k asymptot
% iEig: number of eigenvalues, starts from 1
% lam1: spherical harmonic index 1, starts from 0
% lam2: spherical harmonic index 2, starts from 0
% mu: spherical harmonic index 3, starts from 0
l=iEig-1;
Wi=l*(l+d-2);
Res1=1;
Res2=1;

if l>lam1
    for j=lam1:(l-1)
        Wj=j*(j+d-2);
        ajp1=sqrt((j+1-mu)*(j+1+mu+d-3)/(2*j+d)/(2*j+d-2));
        Res1=Res1*ajp1/(Wi-Wj);
    end
    Res1=Res1*(-1i*K)^(l-lam1);
elseif l<lam1
    for j=(l+1):lam1
        Wj=j*(j+d-2);
        ajp1=sqrt((j-mu)*(j+mu+d-3)/(2*j+d-2)/(2*j+d-4));
        Res1=Res1*ajp1/(Wi-Wj);
    end
    Res1=Res1*(-1i*K)^(lam1-l);
end

if l>lam2
    for j=lam2:(l-1)
        Wj=j*(j+d-2);
        ajp1=sqrt((j+1-mu)*(j+1+mu+d-3)/(2*j+d)/(2*j+d-2));
        Res2=Res2*ajp1/(Wi-Wj);
    end
    Res2=Res2*(-1i*K)^(l-lam2);
elseif l<lam2
    for j=(l+1):lam2
        Wj=j*(j+d-2);
        ajp1=sqrt((j-mu)*(j+mu+d-3)/(2*j+d-2)/(2*j+d-4));
        Res2=Res2*ajp1/(Wi-Wj);
    end
    Res2=Res2*(-1i*K)^(lam2-l);
end

% alternatively use following if lam2=0, mu=0;
% for j=0:(l-1)
%     Wj=j*(j+d-2);
%     ajp1=sqrt((j+1)*(j+1+d-3)/(2*j+d)/(2*j+d-2));
%     Res2=Res2*ajp1^2/(Wi-Wj)^2;
% end
% Res2=Res2*K^(2*l)*(-1)^l;
% Res2=sqrt(Res2);

GLK=Res1*Res2;
end

function GLMK=glm(K,EigK,ORDEig,ORDL,NumLayer,ORDMU)
% output GLM(lamda1+1,lamda2+1,mu+1,l+1)
% this function does not contain the small K limit
% K is the magnitude
% EigK: eigenvalues associated with K.  This is a ORDEig x ORDL matrix.
GLMK=zeros(ORDL,ORDL,ORDMU,ORDEig);
MIN=1e-10;

if abs(K)<MIN
    for M=0:(ORDMU-1)
        for L1=M:(ORDL-1)
            GLMK(L1+1,L1+1,M+1,:)=1;
        end
    end
else
 
for M=0:(ORDMU-1)
    AL=zeros(NumLayer,1);
    for iEig=1:ORDEig
        
        WP=zeros(NumLayer,1);
        WM=zeros(NumLayer,1);
        dJp=zeros(NumLayer,1);
        dJm=zeros(NumLayer,1);
        WPPROD=zeros(NumLayer,1);
        
        LP=NumLayer-1;
        AL(NumLayer)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
        WP(NumLayer)=1/(EigK(iEig,M+1)+LP*(LP+1)); % Q.J.M. changed 8/3/15
        
        LM=M;
        if EigK(iEig,M+1)+LM*(LM+1)==0 % Q.J.M. changed 8/3/15
            WM(M+1)=0;
        else
            WM(M+1)=1/(EigK(iEig,M+1)+LM*(LM+1)); % Q.J.M. changed 8/3/15
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
            
            PL=EigK(iEig,M+1)+LP*(LP+1); % Q.J.M. changed 8/3/15
            AL(IP)=sqrt((LP-M)*(LP+M))/sqrt(4*LP^2-1);
            WP(IP)=1/(PL+WP(IP+1)*(AL(IP+1)*K)^2);
            dJp(IP)=1-dJp(IP+1)*(AL(IP+1)*K*WP(IP+1))^2;
            
            WPPROD(IP)=1i*K*AL(IP)*WP(IP);
            
            PL=EigK(iEig,M+1)+LM*(LM+1); % Q.J.M. changed 8/3/15
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

%{
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
%}

