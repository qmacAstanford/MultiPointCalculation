clear
% close all

chkei=0;
chkres=1;
chkglk=0;
chkglkm=0;

%dimension
d=3.0;

%wavevectors
NK=300;
K0=1e-3;
KF=1e5;
K=transpose(logspace(log10(K0),log10(KF),NK));

ORD=10;                             %number of eigenvalues
ORDL=3;                            %number of spherical harmonic l
ResLayer=500;                     %number of residual layers

R=zeros(NK,ORD);                %eigenvalues
G=zeros(NK,ORD);                %residuals (prefactors)
GLK=zeros(NK,ORD,ORDL);    %residuals (prefactors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate eigenvalues
for I=1:NK
    EI=MatRoots(K(I),d,ORD);
    R(I,:)=transpose(EI);    
end

for I=1:2:ORD
    R(:,I)=real(R(:,I))+1i*abs(imag(R(:,I)));
    R(:,I+1)=real(R(:,I+1))-1i*abs(imag(R(:,I+1)));    
end

if (chkei)
    figure(1);hold
    figure(2);hold
    
    %plot eigenvalues
    for I=1:ORD
        COL=(I-1)/(ORD-1);
        figure(1);plot(K,real(R(:,I)),'-','LineWidth',2,'Color',[COL 0 1-COL])
        figure(2);plot(K,imag(R(:,I)),'-','LineWidth',2,'Color',[COL 0 1-COL])
    end
    figure(1);set(gca,'xscale','log')
    figure(2);set(gca,'xscale','log')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate residuals
if (chkres)
    for I=1:NK
        I
        G(I,:)=transpose(CalRes(R(I,:),K(I),d,ResLayer));
    end

    %plot residuals
    figure(3);hold
    figure(4);hold
    
    for I=1:ORD
        COL=(I-1)/(ORD-1);
        figure(3);plot(K,abs(real(G(:,I))),'-','LineWidth',2,'Color',[COL 0 1-COL])
        figure(4);plot(K,abs(imag(G(:,I))),'-','LineWidth',2,'Color',[COL 0 1-COL])
        %figure(4);plot(K,abs((G(:,I))),'-','LineWidth',2,'Color',[COL 0 1-COL])
    end
    figure(3);set(gca,'xscale','log');set(gca,'yscale','log')
    figure(4);set(gca,'xscale','log');set(gca,'yscale','log')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate glk
if (chkglk)
    for I=1:NK
        I
        GLK(I,:,:)=gl(R(I,:),K(I),d,ORD,ORDL,ResLayer);
    end

    %plot GLKs
    GLKind=2;
    
    figure(5);hold;set(gca,'fontsize',15);leg1=[];
    figure(6);hold;set(gca,'fontsize',15);leg2=[];
    for I=1:ORD
        COL=(I-1)/(ORD-1);

        figure(5);plot(K,real(GLK(:,I,GLKind)),'.-','LineWidth',2,'Color',[COL 0 1-COL])
        figure(6);plot(K,imag(GLK(:,I,GLKind)),'.-','LineWidth',2,'Color',[COL 0 1-COL])
        leg1=[leg1,{['l=',num2str(I)]}];
        leg2=[leg2,{['l=',num2str(I)]}];
    end
    figure(5);set(gca,'xscale','log');set(gca,'yscale','log')
    figure(6);set(gca,'xscale','log');set(gca,'yscale','log')
    
    %using small K expansion
    for I=1:NK
        I
        GLK(I,:,:)=glsmall(R(I,:),K(I),d,ORD,ORDL,ResLayer);
    end

    %plot GLKs
    for I=1:ORD
        COL=(I-1)/(ORD-1);

%         figure(5);plot(K,real(GLK(:,I,GLKind)),'.-','LineWidth',2,'Color',[COL 0 1-COL])
%         figure(6);plot(K,imag(GLK(:,I,GLKind)),'.-','LineWidth',2,'Color',[COL 0 1-COL])
        figure(5);plot(K,real(GLK(:,I,GLKind)),'k.','LineWidth',2)
        figure(6);plot(K,imag(GLK(:,I,GLKind)),'k.','LineWidth',2)
    end
    leg1=[leg1,{'small k asymptot'}];
    leg2=[leg2,{'small k asymptot'}];
    figure(5);set(gca,'xscale','log');set(gca,'yscale','log')
    figure(6);set(gca,'xscale','log');set(gca,'yscale','log')
    figure(5);xlabel('K');ylabel('Real(G_l)');title(['Spherical Harmonic index \lambda=',num2str(GLKind)]);legend(leg1);xlim([K0,KF])
    figure(6);xlabel('K');ylabel('Imag(G_l)');title(['Spherical Harmonic index \lambda=',num2str(GLKind)]);legend(leg2);xlim([K0,KF])
end