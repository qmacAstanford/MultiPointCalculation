function GLK=glsmall(R,k,d,ORDEig,ORDL,ResLayer)

% calculate the residual of all eigenvalues

ResThreshold=1e-10;     % threshold to go from small k asymptot to matrix method
ImagThreshold=1e-8;     % threshold to go from small k asymptot to matrix method
RealThreshold=1e-8;     % threshold to go from small k asymptot to matrix method
NR=length(R);                 % number of eigenvalues (roots)

% get the residues for all roots given in R using recursive relation for
% derivative

GLK=zeros(ORDEig,ORDL); % return value
AL=zeros(ResLayer,1);       % residual layers
lam=0:ORDL-1;                % spherical harmonic L1

for iEig=1:NR                   % eigenvalue index
    for iLam1=1:ORDL        % spherical harmonic index of L1
        lam1=lam(iLam1);
        
        % use asymptotic residual to figure out whether calculate the residual
        % using continued fraction or stay with the asymptotic form for small k
        % limit

        % residual using perturbation expansion for small k
        GLK(iEig,iLam1)=SmallAsympRes(k,iEig,lam1,0,0,d);
%         GLK(iEig,iLam1)=1;
        ResThreshold=1e30;
        
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
        
%         %only keep real part if imaginery part is small
%         if abs(imag(R(iEig)))<ImagThreshold
%             GLK(iEig,iLam1)=real(GLK(iEig,iLam1));
%         end
%         
%         if abs(real(R(iEig)))<RealThreshold
%             GLK(iEig,iLam1)=imag(GLK(iEig,iLam1));
%         end
        
    end
end
end