function GLK=SmallAsympRes(K,iEig,lam1,lam2,mu,d)
% calculate the residual using small k asymptot
% iEig: number of eigenvalues, starts from 1
% lam1: spherical harmonic index 1, starts from 0
% lam2: spherical harmonic index 2, starts from 0
% mu: spherical harmonic index 3, starts from 0
quinnsEdits=1;
if quinnsEdits
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
else
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
end