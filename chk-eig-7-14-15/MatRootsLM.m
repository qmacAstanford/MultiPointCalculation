function Eig=MatRootsLM(k,N)
% output demensions: (l,m)
%

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

