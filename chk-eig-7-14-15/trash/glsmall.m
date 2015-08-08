function GLK=glsmall(R,k,d,ORDEig,ORDL,ResLayer)

% calculate the residual of all eigenvalues

ResThreshold=1e-4;     % threshold to go from small k asymptot to matrix method
ImagThreshold=1e-8;     % threshold to go from small k asymptot to matrix method
NR=length(R);                 % number of eigenvalues (roots)

% get the residues for all roots given in R using recursive relation for
% derivative

GLK=zeros(ORDEig,ORDL); % return value
AL=zeros(ResLayer,1);       % residual layers
lam=0:ORDL-1;               % spherical harmonic L1

for iEig=1:NR                   % eigenvalue index
    for iLam1=1:ORDL        % spherical harmonic index of L1
        lam1=lam(iLam1);
        
        % use asymptotic residual to figure out whether calculate the residual
        % using continued fraction or stay with the asymptotic form for small k
        % limit

        % residual using perturbation expansion for small k
        GLK(iEig,iLam1)=SmallAsympRes(k,iEig,lam1,0,0,d);
        
        if abs(GLK(iEig,iLam1))>ResThreshold

            % residual using continued fraction for large k
            WP=zeros(ResLayer,1);
            dJp=zeros(ResLayer,1);
            WPPROD=zeros(ResLayer,1);

            n=ResLayer-1;
            AL(ResLayer)=n/sqrt(4*n^2-1);
            WP(ResLayer)=1/(R(iEig)+n*(n+1));
            dJp(ResLayer)=1;
            WPPROD(ResLayer)=1i*k*AL(ResLayer)*WP(ResLayer);

            for n=ResLayer-2:-1:0    
                PL=R(iEig)+n*(n+1);
                AL(n+1)=n/sqrt(4*n^2-1);

                WP(n+1)=1/(PL+WP(n+2)*(AL(n+2)*k)^2);
                dJp(n+1)=1-dJp(n+2)*(AL(n+2)*k*WP(n+2))^2;
                WPPROD(n+1)=1i*k*AL(n+1)*WP(n+1);
            end

            if lam1==1
                GLK(iEig,iLam1)=1/dJp(1);  % lam1=0 term
            else
                GLK(iEig,iLam1)=prod(WPPROD(2:(lam1+1)))/dJp(1);  % lam1 >= 1
            end
            
        end
        
    end
end
end

function GLK=SmallAsympRes(K,iEig,lam1,lam2,mu,d)
% calculate the residual using small k asymptot
% iEig: number of eigenvalues
% lam1: spherical harmonic index 1
% lam2: spherical harmonic index 2
% mu: spherical harmonic index 3

l=iEig-1;
mu=mu+1;    % to be consistent with mu=0 case

Wi=l*(l+d-2);
Res1=1;
Res2=1;

if l>lam1
    for j=lam1:(l-1)
        Wj=j*(j+d-2);
        ajp1=sqrt((j+1)*(j+1+d-3)/(2*j+d)/(2*j+d-2));
        Res1=Res1*ajp1^mu/(Wi-Wj);
    end
    Res1=Res1*(-1i*K)^(l-lam1);
elseif l<lam1
    for j=(l+1):lam1
        Wj=j*(j+d-2);
        ajp1=sqrt(j*(j+d-3)/(2*j+d-2)/(2*j+d-4));
        Res1=Res1*ajp1^mu/(Wi-Wj);
    end
    Res1=Res1*(-1i*K)^(lam1-l);
end

if l>lam2
    for j=lam2:(l-1)
        Wj=j*(j+d-2);
        ajp1=sqrt((j+1)*(j+1+d-3)/(2*j+d)/(2*j+d-2));
        Res2=Res2*ajp1^mu/(Wi-Wj);
    end
    Res2=Res2*(-1i*K)^(l-lam2);
elseif l<lam2
    for j=(l+1):lam2
        Wj=j*(j+d-2);
        ajp1=sqrt(j*(j+d-3)/(2*j+d-2)/(2*j+d-4));
        Res2=Res2*ajp1^mu/(Wi-Wj);
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