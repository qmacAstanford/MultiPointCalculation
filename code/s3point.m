function [s3aaa,Qmag,Thet]=s3point
%this function calculates the cubic coefficient
%in free energy expansion

%which chain? (1=Gaussian, 2=WLC, 3=Rigid Rod)
CHAIN=1;
ORDEig=5;
ORDL=5;
NumLayer=500;

%parameters :: polymer
N=10;
NM=10;
FA=0.5;
LAM=-.75;
d=3;

%parameters :: wave vector
NQ=100;NT=100;
Qmag=logspace(-1,2,NQ)/N/NM;
pert=0.;
Thet=linspace(-2*pi+pert,2*pi-pert,NT);
s3aaa=zeros(NQ,NT);
s3aab=zeros(NQ,NT);
s3aba=zeros(NQ,NT);
s3baa=zeros(NQ,NT);

for I=1:NT
    I
    for J=1:NQ
        T=Thet(I);
        Q=Qmag(J);
        
        s3=s3calc(Q,T,N,NM,LAM,FA,d,CHAIN,ORDEig,ORDL,NumLayer);
        s3aaa(I,J)=s3(1,1,1);
        s3aab(I,J)=s3(1,1,2);
        s3aba(I,J)=s3(1,2,1);
        s3baa(I,J)=s3(2,1,1);
    end
end

%make a plot
figure;
if NT==1
    set(gca,'fontsize',15)
    plot(Qmag*N*NM,s3aaa)
    xlabel('Q')
    zlabel('\Gamma')
else
    set(gca,'fontsize',15)
    surf(Qmag*N*NM,Thet,s3aaa,'edgecolor','none','LineStyle','none','FaceLighting','phong');
    set(gca,'xscale','log')
    view([-90,90])
    xlabel('kL')
    ylabel('\theta')
    ylim([min(Thet),max(Thet)])
    zlabel('S_{123}')
    set(gca, 'CLim', [0,0.2]);
    colorbar
end
end

function s3=s3calc(Q,T,N,NM,LAM,FA,d,CHAIN,ORDEig,ORDL,NumLayer)
    %wave vectors
    Q1=Q*[1,0,0];
    Q2=transpose(rotz(T)*Q1(1:3)');
    Q3=-Q1-Q2;
    
    if CHAIN==1
        s3 = s3gaussian(N,NM,LAM,FA,Q1,Q2,Q3,d);
    elseif CHAIN==2
        s3 = s3wlc(N,NM,LAM,FA,Q1,Q2,Q3,ORDEig,ORDL,NumLayer);
    elseif CHAIN==3
        s3 = s3rigid(N,NM,LAM,FA,Q1,Q2,Q3,1);
    end
    s3=s3/power(N*NM,3);
    
end