%Entropy production three dimensional optical trap
%Author Adrian Guel November 2021

clear all;
close all;
clc


fig=figure;
set(fig, 'Position',  [100,100,683,383])
set(gcf,'color','w');

ax1 = subplot(2,2,1);
%set(ax1,'YScale','log')
%set(ax1,'XScale','log')
hold(ax1,'on')
grid(ax1,'on')
colororder(ax1,{'k','b'})
yyaxis(ax1,'left')
ylabel(ax1,'$\mathcal{E}$','Interpreter','Latex','FontSize', 14)
yyaxis(ax1,'right')
ylabel(ax1,'$D_{nn}(t)$','Interpreter','Latex','FontSize', 14)
xlabel(ax1,'$t$','Interpreter','Latex','FontSize', 14)

ax2 = subplot(2,2,2);
%set(ax2,'YScale','log')
%set(ax2,'XScale','log')
hold(ax2,'on')
grid(ax2,'on')
colororder(ax2,{'k','b'})
yyaxis(ax2,'left')
ylabel(ax2,'$\mathcal{E}-\mathcal{E}_m$','Interpreter','Latex','FontSize', 14)
yyaxis(ax2,'right')
ylabel(ax2,'$D_{nn}(t)$','Interpreter','Latex','FontSize', 14)
xlabel(ax2,'$t$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(2,2,3);
set(ax3,'YScale','log')
%set(ax3,'XScale','log')
hold(ax3,'on')
grid(ax3,'on')
colororder(ax3,{'k','b'})
yyaxis(ax3,'left')
ylabel(ax3,'$\Pi_i$','Interpreter','Latex','FontSize', 14)
yyaxis(ax3,'right')
ylabel(ax3,'$D_{nn}(t)$','Interpreter','Latex','FontSize', 14)
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)

ax4 = subplot(2,2,4);
%set(ax4,'YScale','log')
%set(ax4,'XScale','log')
hold(ax4,'on')
grid(ax4,'on')
colororder(ax4,{'k','b'})
yyaxis(ax4,'left')
ylabel(ax4,'$\Pi-\Pi_i$','Interpreter','Latex','FontSize', 14)
yyaxis(ax4,'right')
ylabel(ax4,'$D_{nn}(t)$','Interpreter','Latex','FontSize', 14)
xlabel(ax4,'$t$','Interpreter','Latex','FontSize', 14)

a=0.1;
b=[0,0,0,1];
tx=[6,6,6,6];
tspan = [0 10];   
%t=0:0.005:2;
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
%tspan=logspace(0,5);
maxn=4;
crls=rand([maxn,3]);

rng('default');
rng(5);
legaux={};
for n=4:maxn
    %A=[0,1;-1,-2];
    %A=[0,0.5,0.5; 0,0,0.5; 0,1,0];
    %positive definite random matrix
    % il=rand(n);
    % A=il*il.';
    %positive definite hurwitz matrix
    sys=rss(n);
    csys = canon(sys,'companion');
    A=csys.A';

    B=flip(csys.B); ninputs=size(B);
    %u=ones(ninputs(2),1);
    u=zeros(ninputs(2),1);

    Sigma0=0.1*eye(n);
    %X0=-5 +10*rand(1,n);
    X0=zeros(1,n); x0(5)=1;
    At=A';
    D=0.01*eye(n);
    y0=Sigma0(:)';
    [t,y] = ode45(@(t,y) Sigma(t,y,A,D,n,tx,b,a), tspan, y0, opts);
    [t,M] = ode45(@(t,y) Mu(t,y,A,B,u,a,tx), t, X0, opts);
    %M=Mu(t,X0,A,n);
    %y=Sigma(t,Sigma0,A,D,n);
    % 
    dSdt=zeros(length(t),1);
    Pi=zeros(length(t),1);
    Phi=zeros(length(t),1);
    E=zeros(length(t),1);
    Pit=zeros(length(t),1);
    Ip=zeros(length(t),1);
    Ep=zeros(length(t),1);
    
    dSkdt=zeros(length(t),n);
    Pik=zeros(length(t),n);
    Eaux=zeros(length(t),n);
    Piksum=zeros(length(t),1);
    DSkdtsum=zeros(length(t),1);
    Esumi=zeros(length(t),1);

    for k=1:length(dSdt)
        Sx=reshape(y(k,:),n,n);
        dSdt(k)=0.5*trace(Sx\(A*Sx+Sx*A'+2*D));
        Pi(k)=((A*M(k,:)')'/D)*(A*M(k,:)')+trace(At/D*A*Sx+Sx\D+2*A);
        %Pit(k)=trace((A*M(k,:))*(A*M(k,:))')*trace(inv(D))+trace(inv(Sx))*trace(Sx*A*Sx*A'+2*Sx*A*D+D^2)*trace(inv(D));
        Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+(1/4)*trace(inv(Sx))*((trace(A*Sx+Sx*At+2*D))^2)*trace(inv(D));
        Phi(k)=((A*M(k,:)')'/D)*(A*M(k,:)')+trace(At/D*A*Sx+A);
        E(k)=(A*M(k,:)')'/Sx*(A*M(k,:)')+0.5*trace((Sx\(A*Sx+Sx*A'+2*D))^2);
        Ip(k)=(1/n)*(trace(inv(Sx))*Pit(k)*trace(D))+dSdt(k)^2;
        Ep(k)=(1/n)*(trace(inv(Sx))*Pi(k)*trace(D))+dSdt(k)^2-2*S2(Sx\D+A);
    end
    
    
    %y(:,[2 3 4 5 6 8 9 10 11 12 14 15 16 17 18 20 21 22 23 24])=[]; %%5order
    y(:,[2 3 4 5 7 8 9 10 12 13 14 15])=[];
    for i=1:n
        dSkdt(:,i)=0.25*((2*A(i,i)*(y(:,i))+2*D(i,i))./(y(:,i))).^2;
        Pik(:,i)=((A(i,i)*M(:,i)).^2)./D(i,i)+(1/4)*((2*A(i,i)*(y(:,i))+2*D(i,i)).^2)./((y(:,i))*D(i,i));
        Eaux(:,i)=((A(i,i)*M(:,i)).^2)./(y(:,i))+0.5*((2*A(i,i)*(y(:,i))+2*D(i,i))./(y(:,i))).^2;
        Piksum=Piksum+Pik(:,i);%*D(i,i))./(y(:,i));
        DSkdtsum=DSkdtsum+dSkdt(:,i);
        Esumi=Esumi+Eaux(:,i);
        legaux{i}=sprintf('i=%d', i);
        legaux2{i}=sprintf('i=%d', i);
    end
    Esumd=DSkdtsum+Piksum;
    
  
    yyaxis(ax1,'left')
    plot(ax1,t,Eaux(:,1),'.-');
    plot(ax1,t,Eaux(:,2),':');
    plot(ax1,t,Eaux(:,3),'--');
    plot(ax1,t,Eaux(:,4),'-');
    plot(ax1,t,E,'r:');
    yyaxis(ax1,'right')
    plot(ax1,t,b(n)*(1/(abs(a)*sqrt(pi))).*exp(-((t-tx(1))./a).^2));
    
    %plot(ax2,t,Ep-E,'Color',[crls(n,:)])
    yyaxis(ax2,'left')    
    plot(ax2,t,E-Esumd);
    yyaxis(ax2,'right')
    plot(ax2,t,b(n)*(1/(abs(a)*sqrt(pi))).*exp(-((t-tx(1))./a).^2));    
    %%5thorder
    %plot(ax3,t,Pik(:,1),t,Pik(:,2),t,Pik(:,3),t,Pik(:,4),t,Pik(:,5),t,Piksum,'-.');
    %plot(ax4,t,Pi-Pik(:,1),t,Pi-Pik(:,2),t,Pi-Pik(:,3),t,Pi-Pik(:,4),t,Pi-Pik(:,5),t,Pi-Piksum,'-.');
    
    %%4thorder
    yyaxis(ax3,'left')
    plot(ax3,t,Pik(:,1),t,Pik(:,2),t,Pik(:,3),t,Pik(:,4),t,Piksum);
    yyaxis(ax3,'right')
    plot(ax3,t,b(n)*(1/(abs(a)*sqrt(pi))).*exp(-((t-tx(1))./a).^2));
    
    yyaxis(ax4,'left')
    plot(ax4,t,Pi-Pik(:,1),t,Pi-Pik(:,2),t,Pi-Pik(:,3),t,Pi-Pik(:,4),t,Pi-Piksum);
    yyaxis(ax4,'right')
    plot(ax4,t,b(n)*(1/(abs(a)*sqrt(pi))).*exp(-((t-tx(1))./a).^2));
    
    legaux{5}='i=m';
    legaux2{5}='i=m';
end

%sprintf('Some text %g', someValue)
leg1 = legend(ax3,legaux);
leg2 = legend(ax4,legaux2);
%set(leg1);

function dydt = Sigma(t,y,A,D,n,tx,b,a)
   At=A';
   aux2=reshape(y,n,n);
   Dx=Dfunction4(D,t,b,a,tx);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*Dx(:);
end
 
function dydt = Mu(t,y,A,B,u,a,tx)
   dydt=A*y+B*(u*(1/(abs(a)*sqrt(pi))).*exp(-((t-tx(1))./a).^2));
end

% function dydt = Sigma(t,y,A,D,n)
%    dydt=zeros(n*n,length(t));
%    P2=zeros(n*n,length(t));
%    for k=1:length(t)
%     p1=expm(A*t(k))*y*expm(A*t(k))';
%         for l=1:k
%             aux=2*expm(A*t(l))*D*expm(A*t(l));
%             P2(:,l)=aux(:)';
%         end
%     dydt(:,k)=p1(:)+ trapz(t(1:k),P2(:,1:k),2);
%    end
% end

% function dydt = Mu(t,y,A,n)
%     dydt=zeros(n,length(t));
%     for k=1:length(t)
%      dydt(:,k)=expm(A*t(k))*y';
%     end
% end
