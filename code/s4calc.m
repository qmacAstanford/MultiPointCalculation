function s4=s4calc(Q,T,N,NM,LAM,FA,d,CHAIN,ORDEig,ORDL,NumLayer)
    %wave vectors
    Q1=Q*[1,0,0];
    Q2=transpose(rotz(T)*Q1(1:3)');
    Q3=-Q1;
    Q4=-Q2;
    
    if CHAIN==1
        s4 = s4gaussian(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d);
    elseif CHAIN==2
        s4 = s4wlc(N,NM,LAM,FA,Q1,Q2,Q3,Q4,ORDEig,ORDL,NumLayer);
    elseif CHAIN==3
        s4 = s4rigid(N,NM,LAM,FA,Q1,Q2,Q3,Q4,1);
    end
    s4=s4/power(N*NM,4);
end