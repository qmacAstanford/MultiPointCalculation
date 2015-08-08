function [val,SepRes]=gkwlc(N,k,d,ORDmax,ORD,ResLayer)

% Calculate the Fourier transform of the Green function
% for the wormlike chain in d-dimensions
%
% Andrew Spakowitz (4/14/15)

% Fill in unset optional values

switch nargin
    case 2
        d=3;
        ORDmax=20;
        ORD=20;
        ResLayer=500;
    case 3
        ORDmax=20;
        ORD=20;
        ResLayer=500;        
    case 4
        ORD=20;
        ResLayer=500;        
    case 5
        ResLayer=500;        
end

% If dimensions is 2, reset value to small perturbation above 2

if d==2
    d=2+1e-10;
end

% Reset N to a row vector if entered as a column

if iscolumn(N)==1
    N=transpose(N);
end

val=zeros(length(k),length(N));
SepRes=zeros(length(k),ORDmax);

% calculate the roots or eigenvalues of the Schrodinger equation
% k is a vector of all frequencies, for each k, get the roots

for j=1:length(k)
    
    % calculate the eigenvalues
    R=MatRoots(k(j),d,ORD);
    NR=ORDmax;

    % get the residues for all roots of each k(j)
    Residue=CalRes(R(1:NR),k(j),d,ResLayer);

    for I=1:NR
        val(j,:)=val(j,:)+exp(R(I)*N)*Residue(I);
        SepRes(j,I)=exp(R(I)*min(N))*Residue(I);
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

function value=Epsilon(k,d,alpha)

% eigenvalues using large k asymptotic expansion
% generates epsilon^{s}_r (see the paper)

I=complex(0,1);
beta=-sqrt(2)/4*(1+I);
m=(d-3)/2;
n=2*k+m+1;

epsilon_0=(-1/2/beta)^(-1)*(n/2);
epsilon_1=(-1/2/beta)^( 0)*(-1/8*(n^2+3-3*m^2)-m*(m+1));
epsilon_2=(-1/2/beta)^( 1)*(-1/2^5*n*(n^2+3-9*m^2));
epsilon_3=(-1/2/beta)^( 2)*(-1/2^8*(5*n^4+34*n^2+9)-(102*n^2+42)*m^2+33*m^4);
epsilon_4=(-1/2/beta)^( 3)*(-1/2^11*n*(33*n^4+410*n^2+405)-(1230*n^2+1722)*m^2+813*m^4);
epsilon_5=(-1/2/beta)^( 4)*(-1/2^12*9*(7*n^6+140*n^4+327*n^2+54-(420*n^4+1350*n^2+286)*m^2+(495*n^2+314)*m^4-82*m^6));

value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha^2+...
      epsilon_4*alpha^3+epsilon_5*alpha^4;

function Res=CalRes(R,k,d,ResLayer)

% calculate the residual of all eigenvalues

ResThreshold=1e-12;     % threshold to go from small k asymptot to matrix method
NR=length(R);           % number of eigenvalues (roots)

% get the residues for all roots given in R using recursive relation for
% derivative
Res=zeros(NR,1);    % residual vector, each row corresponds to each eigenvalue

% find the residual
for i=1:NR

    % use asymptotic residual to figure out whether calculate the residual
    % using continued fraction or stay with the asymptotic form for small k
    % limit
    Res(i)=SmallAsympRes(k,i,d);
    if abs(Res(i))>ResThreshold
        if k<=1

            % residual using continued fraction for small k
            p=R(i);
            W=p+(ResLayer+d-2)*ResLayer;
            Wprime=1;
            for L=ResLayer:-1:1
                AL=k*sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));         % d-dimensional case
                Wprime=1-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm+AL^2/W;
            end
            Res(i)=1/Wprime;
        else

            % residual using continued fraction for large k
            p=R(i);
            W=(p+(ResLayer+d-2)*ResLayer)/k;
            Wprime=1/k;
            for L=ResLayer:-1:1
                AL=sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));           % d-dimensional case
                Wprime=1/k-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm/k+AL^2/W;
            end
            Res(i)=1/(k*Wprime);
        end
    end
end

function Res=SmallAsympRes(K,n,d)

% calculate the residual using small k asymptot
n=n-1;
Res=1;
Wn=-n*(n+1);

for j=1:n
    Wjm1=-(j-1)*j;
    aj=sqrt(j*(j+d-3)/(2*j+d-2)/(2*j+d-4));
    Res=Res*aj/(Wn-Wjm1);
end

Res=Res^2*K^(2*n)*(-1)^n;