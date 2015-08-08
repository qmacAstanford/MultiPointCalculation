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
if k<=1
    Eig(:,M+1)=-TempMat(1:N);
else
    Eig(:,M+1)=-TempMat(1:N)*k;
end

end


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

end