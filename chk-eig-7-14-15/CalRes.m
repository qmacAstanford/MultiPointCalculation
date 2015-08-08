function Res=CalRes(R,k,d,ResLayer)

% calculate the residual of all eigenvalues

ResThreshold=1e-12;     % threshold to go from small k asymptot to matrix method
ImagThreshold=1e-8;     % threshold to go from small k asymptot to matrix method
NR=length(R);           % number of eigenvalues (roots)

% get the residues for all roots given in R using recursive relation for
% derivative
Res=zeros(NR,1);    % residual vector, each row corresponds to each eigenvalue

% find the residual
for n=1:NR

    % use asymptotic residual to figure out whether calculate the residual
    % using continued fraction or stay with the asymptotic form for small k
    % limit

    Res(n)=SmallAsympRes(k,n,d);
    if abs(Res(n))>ResThreshold
        if k<=1

            % residual using continued fraction for small k
            p=R(n);
            W=p+(ResLayer+d-2)*ResLayer;
            Wprime=1;
            for L=ResLayer:-1:1
                AL=k*sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));         % d-dimensional case
                Wprime=1-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm+AL^2/W;
            end
            Res(n)=1/Wprime;
        else

            % residual using continued fraction for large k
            p=R(n);
            W=(p+(ResLayer+d-2)*ResLayer)/k;
            Wprime=1/k;
            for L=ResLayer:-1:1
                AL=sqrt(L*(L+d-3)/(2*L+d-2)/(2*L+d-4));           % d-dimensional case
                Wprime=1/k-AL^2*Wprime/W^2;
                PLm=p+(L+d-2-1)*(L-1);                            % d-dimensional case
                W=PLm/k+AL^2/W;
            end
            Res(n)=1/(k*Wprime);
        end
    end
    if abs(imag(R(n)))<ImagThreshold
        Res(n)=real(Res(n));
    end
end

function Res=SmallAsympRes(K,n,d)

% calculate the residual using small k asymptot

l=n-1;
Res=1;
Wl=l*(l+d-2);

for j=0:(l-1)
    Wj=j*(j+d-2);
    ajp1=sqrt((j+1)*(j+1+d-3)/(2*j+d)/(2*j+d-2));
    Res=Res*ajp1^2/(Wl-Wj)^2;
end

Res=Res*K^(2*l)*(-1)^l;