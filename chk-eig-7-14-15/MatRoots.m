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