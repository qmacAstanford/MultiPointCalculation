function Eig=LargeKEigenValues(k,ORDEig,ORDL,d)
% Returs an ORDEig x ORDEig  where entries are for l+1,mu+1
% where l is the order of the eigenvalue
% and mu is is the azimuthal quantum number (i.e. m)

Eig=zeros(ORDEig,ORDL);
mu=0:(ORDL-1);

for I=1:floor(ORDEig/2)
    l=I-1;
    alpha=1/sqrt(8*k);
    Eig(2*l+1,:)=1i*k-mu.*(mu+d-2)-Epsilon(l,d,alpha,mu); % Q.J.M. changes this line 8/1/15
    Eig(2*l+2,:)=conj(Eig(2*l+1,:));
end

    
end

function value=Epsilon(l,d,alpha,mu)

% eigenvalues using large k asymptotic expansion
% generates epsilon^{s}_r (see the paper)

I=complex(0,1);
beta=-sqrt(2)/4*(1+I);
m=mu+(d-3)/2;  % Q.J.M. changes this line 8/1/15
n=2*l+m+1;  % also know as s

epsilon_0=(-1/2/beta)^(-1)*(n/2);
epsilon_1=(-1/2/beta)^( 0)*(-1/8*(n.^2+3-3*m.^2)-m.*(m+1));
epsilon_2=(-1/2/beta)^( 1)*(-1/2^5*n.*(n.^2+3-9*m.^2));
epsilon_3=(-1/2/beta)^( 2)*(-1/2^8*(5*n.^4+34*n.^2+9)-(102*n.^2+42).*m.^2+33*m.^4);
epsilon_4=(-1/2/beta)^( 3)*(-1/2^11*n.*(33*n.^4+410*n.^2+405)-(1230*n.^2+1722).*m.^2+813*m.^4);
epsilon_5=(-1/2/beta)^( 4)*(-1/2^12*9*(7*n.^6+140*n.^4+327*n.^2+54-(420*n.^4+1350*n.^2+286).*m.^2+(495*n.^2+314).*m.^4-82*m.^6));

value=epsilon_0/alpha+epsilon_1+epsilon_2*alpha+epsilon_3*alpha^2+...
      epsilon_4*alpha^3+epsilon_5*alpha^4;
end