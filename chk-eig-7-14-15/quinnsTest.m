
clear
clc
testMatRoots=0;
chooseCutoff=0;
plotEig=0;
prefactor=1;

if prefactor
    manypoints=1;
    if manypoints
        npts=100;
    else
        npts=1;
    end
    ORDEig=4;
    ORDL=5;

    R=zeros(npts,ORDEig,ORDL); % which k, which eigenvalue, which m
    if manypoints
        K=logspace(-3,5,npts);
    else
        K=50;
    end
    
    % R( K index, eig index,m index )
    for count=1:npts
        if K(count) < 200
            R(count,:,:)=IntermediateKEigenValues(K(count),ORDEig,ORDL);
        else
            R(count,:,:)=LargeKEigenValues(K(count),ORDEig,ORDL,3);
        end
    end
    
%     figure(10)
%     hold on
%     semilogx(K,R(:,1,1),'o',K,R(:,2,1),'o',K,R(:,3,1),'o')
    GLMK_largeKV2=zeros(npts,ORDL,ORDL,ORDL,ORDEig)*NaN;
    GLMK_largeK=zeros(npts,ORDL,ORDL,ORDL,ORDEig)*NaN;
    GLMK_smallK=zeros(npts,ORDL,ORDL,ORDL,ORDEig)*NaN;
    usinggl=zeros(npts,ORDL,ORDEig);
    usingCalRes=zeros(npts,ORDEig);
    for count=1:npts
        GLMK_largeKV2(count,:,:,:,:)=glmV2(K(count),squeeze(R(count,:,:)),ORDEig,ORDL,500);  
        GLMK_largeK(count,:,:,:,:)=  glm(K(count),squeeze(R(count,:,:)),ORDEig,ORDL,501);    
        
        for lam1=0:ORDL-1
            for lam2=0:ORDL-1
                for iEig=1:ORDEig
                    for mu=0:min(lam1,lam2)
                        GLMK_smallK(count,lam1+1,lam2+1,mu+1,iEig)=SmallAsympRes(K(count),iEig,lam1,lam2,mu,3);
                    end
                end
            end
        end
        usinggl(count,:,:)=gl(R(count,:,1),K(count),3,ORDEig,ORDL,500)';%  zeros(ORDEig,ORDL);
        usingCalRes(count,:)=CalRes(R(count,:,1),K(count),3,500);
    end
    doesItWorkForLowL=0;
    if doesItWorkForLowL
        figure(1)

        lam1=0;
        lam2=0;
        mu=0;
        L=0;
        hold on
        GLMK_largeK(:,lam1+1,lam2+1,mu+1,L+1);
        semilogx(K,real(GLMK_largeK(:,lam1+1,lam2+1,mu+1,L+1)),'-x',...
                 K,real(GLMK_smallK(:,lam1+1,lam2+1,mu+1,L+1)),...
                 K,real(    usinggl(:,lam1+1,            L+1)),'.',...
                 K,real(usingCalRes(:,                   L+1)),'o');
        legend('largK','smallk','old')
        title(sprintf('Compair methodes lam1=%d, lam2=%d, mu=%d, L=%d',lam1,lam2,mu,L))
        xlabel('K')
        ylim([-2,2])
        hold off
        
        figure(2)
        hold on
          loglog(K,abs(real(GLMK_largeK(:,lam1+1,lam2+1,mu+1,L+1))),...
                 K,abs(real(GLMK_smallK(:,lam1+1,lam2+1,mu+1,L+1))),...
                 K,abs(real(    usinggl(:,lam1+1,            L+1))),'.',...
                 K,abs(real(usingCalRes(:,                   L+1))),'o');
        legend('largK','smallk','old')
        title(sprintf('Compair methodes lam1=%d, lam2=%d, mu=%d, L=%d',lam1,lam2,mu,L))
        xlabel('K')  
        hold off  
        figure(1);set(gca,'xscale','log');
        figure(2);set(gca,'xscale','log');set(gca,'yscale','log')   
    else    
        figure(1)

        lam1=3;
        lam2=3;
        mu=0;
        L=2;
        hold on
        semilogx(K,real(GLMK_largeK(:,lam1+1,lam2+1,mu+1,L+1)),...
                 K,real(GLMK_largeKV2(:,lam1+1,lam2+1,mu+1,L+1)),'x',...
                 K,real(GLMK_smallK(:,lam1+1,lam2+1,mu+1,L+1)));
        legend('largK','v2','smallk')
        title(sprintf('Re(G_l_m(K)) lam1=%d, lam2=%d, mu=%d, L=%d',lam1,lam2,mu,L))
        xlabel('K')
        ylim([-2,2])
        hold off

        figure(2)
        hold on
          loglog(K,abs(real(GLMK_largeK(:,lam1+1,lam2+1,mu+1,L+1))),...
                 K,abs(real(GLMK_largeKV2(:,lam1+1,lam2+1,mu+1,L+1))),'x',...
                 K,abs(real(GLMK_smallK(:,lam1+1,lam2+1,mu+1,L+1))));
        legend('largK','v2','smallk')
        title(sprintf('Compair methodes lam1=%d, lam2=%d, mu=%d, L=%d',lam1,lam2,mu,L))
        xlabel('K')  
        hold off  
        figure(1);set(gca,'xscale','log');
        figure(2);set(gca,'xscale','log');set(gca,'yscale','log')     
    end 
    
    
    
    
end

if plotEig
    npts=70;
    ORDEig=4;
    ORDL=5;

    R=zeros(npts,ORDEig,ORDL);
    K=logspace(-1,5,npts);

    for count=1:npts
        if K(count) < 200
            R(count,:,:)=IntermediateKEigenValues(K(count),ORDEig,ORDL);
        else
            R(count,:,:)=LargeKEigenValues(K(count),ORDEig,ORDL,3);
        end
    end

    figure
    semilogx(K,real(R(:,1,1)),'r',...
             K,real(R(:,1,2)),'r.',...
             K,real(R(:,1,3)),'ro',...
             K,real(R(:,1,4)),'r--',...
             K,real(R(:,2,1)),'b',...
             K,real(R(:,2,2)),'b.',...
             K,real(R(:,2,3)),'bo',...
             K,real(R(:,2,4)),'b--',...         
             K,real(R(:,3,1)),'g',...
             K,real(R(:,3,2)),'g.',...
             K,real(R(:,3,3)),'go',...
             K,real(R(:,3,4)),'g--',...         
             K,real(R(:,4,1)),'k',...
             K,real(R(:,4,2)),'k.',...
             K,real(R(:,4,3)),'ko',...
             K,real(R(:,4,4)),'k--');
    xlabel('|K|')
    ylabel('Re(Eigenvalue)')
    legend('location','eastoutside',...
           'L=0, M=0','L=0, M=1','L=0, M=2','L=0, M=3',...
           'L=1, M=0','L=1, M=1','L=1, M=2','L=1, M=3',...
           'L=2, M=0','L=2, M=1','L=2, M=2','L=2, M=3',...
           'L=3, M=0','L=3, M=1','L=3, M=2','L=3, M=3') 
    title('Cutoff is at K=200')

    % Now imagionary
    figure
    loglog(...
             K,imag(R(:,1,1)),'r',...
             K,imag(R(:,1,2)),'r.',...
             K,imag(R(:,1,3)),'ro',...
             K,imag(R(:,1,4)),'r--',...         
             K,imag(R(:,3,1)),'g',...
             K,imag(R(:,3,2)),'g.',...
             K,imag(R(:,3,3)),'go',...
             K,imag(R(:,3,4)),'g--'...                 
             );   
    ylabel('Im(Eigenvalue)')
    legend('location','eastoutside',...
           'L=0, M=0 intr','L=0, M=1 intr','L=0, M=2 intr','L=0, M=3 intr',...
           'L=2, M=0 intr','L=2, M=1 intr','L=2, M=2 intr','L=2, M=3 intr'...
           ); 
end
if chooseCutoff
    npts=70;
    ORDEig=4;
    ORDL=5;

    R=zeros(npts,ORDEig,ORDL);
    K=logspace(-1,5,npts);

    % now intermediate
    for count=1:npts
        R(count,:,:)=IntermediateKEigenValues(K(count),ORDEig,ORDL);
    end
    Rinter=R;
    figure
    semilogx(K,real(R(:,1,1)),'r',...
             K,real(R(:,1,2)),'r.',...
             K,real(R(:,1,3)),'ro',...
             K,real(R(:,1,4)),'r--',...
             K,real(R(:,2,1)),'b',...
             K,real(R(:,2,2)),'b.',...
             K,real(R(:,2,3)),'bo',...
             K,real(R(:,2,4)),'b--',...         
             K,real(R(:,3,1)),'g',...
             K,real(R(:,3,2)),'g.',...
             K,real(R(:,3,3)),'go',...
             K,real(R(:,3,4)),'g--',...         
             K,real(R(:,4,1)),'k',...
             K,real(R(:,4,2)),'k.',...
             K,real(R(:,4,3)),'ko',...
             K,real(R(:,4,4)),'k--');
    xlabel('|K|')
    ylabel('Re(Eigenvalue)')
    legend('location','eastoutside',...
           'L=0, M=0','L=0, M=1','L=0, M=2','L=0, M=3',...
           'L=1, M=0','L=1, M=1','L=1, M=2','L=1, M=3',...
           'L=2, M=0','L=2, M=1','L=2, M=2','L=2, M=3',...
           'L=3, M=0','L=3, M=1','L=3, M=2','L=3, M=3') 

    % Now large K  
    R=zeros(npts,ORDEig,ORDL);
    for count=1:npts
        R(count,:,:)=LargeKEigenValues(K(count),ORDEig,ORDL,3);
    end
    Rlong=R;

    hold on
    semilogx(K,real(R(:,1,1)),'rx',...
             K,real(R(:,1,2)),'r-x',...
             K,real(R(:,1,3)),'r.-',...
             K,real(R(:,1,4)),'rx',...        
             K,real(R(:,3,1)),'gx',...
             K,real(R(:,3,2)),'g-x',...
             K,real(R(:,3,3)),'g.-',...
             K,real(R(:,3,4)),'gx');
    xlabel('|K|')
    ylabel('Re(Eigenvalue)')
    title('Real comparison')



    % legend('location','eastoutside',...
    %        'L=0, M=0','L=0, M=1','L=0, M=2','L=0, M=3',...
    %        'L=1, M=0','L=1, M=1','L=1, M=2','L=1, M=3',...
    %        'L=2, M=0','L=2, M=1','L=2, M=2','L=2, M=3',...
    %        'L=3, M=0','L=3, M=1','L=3, M=2','L=3, M=3')

     ylim([-150,0])

     %  --- figure out which is better where ----
    figure 
    R=abs(real(R-Rinter));

    loglog(K,real(R(:,1,1)),'r',...
             K,real(R(:,1,2)),'r.',...
             K,real(R(:,1,3)),'ro',...
             K,real(R(:,1,4)),'r--',...
             K,real(R(:,2,1)),'b',...
             K,real(R(:,2,2)),'b.',...
             K,real(R(:,2,3)),'bo',...
             K,real(R(:,2,4)),'b--',...         
             K,real(R(:,3,1)),'g',...
             K,real(R(:,3,2)),'g.',...
             K,real(R(:,3,3)),'go',...
             K,real(R(:,3,4)),'g--',...         
             K,real(R(:,4,1)),'k',...
             K,real(R(:,4,2)),'k.',...
             K,real(R(:,4,3)),'ko',...
             K,real(R(:,4,4)),'k--');
    xlabel('|K|')
    ylabel('abs(Re(intermediateK)-Re(largeK))')
    legend('location','eastoutside',...
           'L=0, M=0','L=0, M=1','L=0, M=2','L=0, M=3',...
           'L=1, M=0','L=1, M=1','L=1, M=2','L=1, M=3',...
           'L=2, M=0','L=2, M=1','L=2, M=2','L=2, M=3',...
           'L=3, M=0','L=3, M=1','L=3, M=2','L=3, M=3')
    title('Choosing real cutoff')

    % Now imagionary
    figure
    loglog(...
             K,imag(Rinter(:,1,1)),'r',...
             K,imag(Rinter(:,1,2)),'r.',...
             K,imag(Rinter(:,1,3)),'ro',...
             K,imag(Rinter(:,1,4)),'r--',...         
             K,imag(Rinter(:,3,1)),'g',...
             K,imag(Rinter(:,3,2)),'g.',...
             K,imag(Rinter(:,3,3)),'go',...
             K,imag(Rinter(:,3,4)),'g--',...         
             K,imag(Rlong(:,1,1)),'b',...
             K,imag(Rlong(:,1,2)),'b.',...
             K,imag(Rlong(:,1,3)),'bo',...
             K,imag(Rlong(:,1,4)),'b--',...         
             K,imag(Rlong(:,3,1)),'k',...
             K,imag(Rlong(:,3,2)),'k.',...
             K,imag(Rlong(:,3,3)),'ko',...
             K,imag(Rlong(:,3,4)),'k--'...         
             );   
    ylabel('Im(Eigenvalue)')
    legend('location','eastoutside',...
           'L=0, M=0 intr','L=0, M=1 intr','L=0, M=2 intr','L=0, M=3 intr',...
           'L=2, M=0 intr','L=2, M=1 intr','L=2, M=2 intr','L=2, M=3 intr',...
           'L=0, M=0 long','L=0, M=1 long','L=0, M=2 long','L=0, M=3 long',...
           'L=2, M=0 long','L=2, M=1 long','L=2, M=2 long','L=2, M=3 long'...
           ); 
    title('imaginary comparison')

    figure()
    R=abs(imag(Rinter-Rlong));
      loglog(K,R(:,1,1),'r',...
             K,R(:,1,2),'r.',...
             K,R(:,1,3),'ro',...
             K,R(:,1,4),'r--',...
             K,R(:,3,1),'b',...
             K,R(:,3,2),'b.',...
             K,R(:,3,3),'bo',...
             K,R(:,3,4),'b--');
    xlabel('|K|')
    ylabel('abs(Im(intermediateK)-Im(largeK))')
    legend('location','eastoutside',...
           'L=0, M=0','L=0, M=1','L=0, M=2','L=0, M=3',...
           'L=2, M=0','L=2, M=1','L=2, M=2','L=2, M=3')
    title('Choosing imaginary cutoff')
end
   
if testMatRoots

    ORDEig=6; %Needs to be even

    npts=100;
    R=zeros(npts,ORDEig);
    K=logspace(-1,2,npts);

    for count=1:npts
        R(count,:)=MatRoots(K(count),3,ORDEig)';
    end


    for I=1:2:ORDEig
        R(:,I)=real(R(:,I))+1i*abs(imag(R(:,I)));
        R(:,I+1)=real(R(:,I+1))-1i*abs(imag(R(:,I+1)));    
    end

    R=R';

    figure
    semilogx(K,real(R(1,:)),'.-',...
             K,real(R(2,:)),'.-',...
             K,real(R(3,:)),...
             K,real(R(4,:)),...
             K,real(R(5,:))...
             )
     xlabel('K')
     ylabel('Output of Matroots')
     legend('l=1','l=2','l=3','l=4','l=5')
     title('reak part, same as figure 1A')

     figure
     semilogx(K,imag(R(1,:)),'.-',...
             K,imag(R(2,:)),'.-',...
             K,imag(R(3,:)),...
             K,imag(R(4,:)),...
             K,imag(R(5,:)),...
             K,imag(R(5,:))...
             )
     xlabel('K')
     ylabel('Output of Matroots')
     legend('l=0','l=1','l=2','l=3','l=4','l=5')
     title('imagionary part, same as figure 1A')
end