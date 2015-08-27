function Eig=SmallEigenvales(k,ORDEig,ORDL,D)

Eig=zeros(ORDEig,ORDL)*NaN;

for l=0:(ORDEig-1)
    for M=0:(ORDL-1)
        b0=-l*(l+D-2);
        b2=alm(l,M,D)^2/(2*l+D-3) - (alm(l+1,M,D)^2)/(2*l+D-1);
        b4=alm(l,M,D)^4/((2*l+D-3)^3) ...
           -alm(l+1,M,D)^4/((2*l+D-1)^3) ...
           +alm(l+1,M,D)^2 * alm(l+2,M,D)^2 /(2*(2*l+D)*(2*l+D-1)^2) ...
           -alm(l,M,D)^2 * alm(l-1,M,D)^2 /(2*(2*l+D-4)*(2*l+D-3)^2) ...
           -2*alm(l,M,D)^2 * alm(l+1,M,D)^2 / ((2*l+D-1)^2 *(2*l+D-3)^2);
       Eig(l+1,M+1)=b0+b2*k^2+b4*k^4;
    end
end

end
function out=alm(l,mu,D)
% you may use a vector for l in which case the output will be the same size
    out=sqrt((l-mu ).*(l+mu+D-3)./((2*l+D-2).*(2*l+D-4)));
    if l(1)==-1
        out(1)=0;
    end
end