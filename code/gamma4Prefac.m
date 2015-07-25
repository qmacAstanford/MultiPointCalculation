function out=gamma4Prefac(Ivec,AlphaSpace,order,FA,n)
% Ivec is binary row vector of length n
% AlphaVec is binary row vector of length n
% order is an inter between 1 and n!

aPath = permute(AlphaSpace,order);
F=[FA,1-FA];
out=F(aPath(1)+1); % f_alpha_1

for binomialNumber=1:n-1
    a_current=aPath(binomialNumber);
    a_next=aPath(binomialNumber+1);
    if Ivec(binomialNumber)==1
        out=out*Delta2(a_current,a_next)*(1-F(a_current+1)); % Delta*(1-f_alpha_1)
    else
        out=out*F(a_next+1); % f_alpha_2
    end
end

end
function out = Delta2(alpha1,alpha2)
    if alpha1==alpha2
        out=1;
    else
        out=-1;
    end
end