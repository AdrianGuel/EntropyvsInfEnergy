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
ylabel(ax1,'$\mathcal{E}$','Interpreter','Latex','FontSize', 14)
xlabel(ax1,'$\dot{S}$','Interpreter','Latex','FontSize', 14)


ax2 = subplot(2,2,[2 4]);
%set(ax2,'YScale','log')
%set(ax2,'XScale','log')
hold(ax2,'on')
grid(ax2,'on')
axis(ax2,'square')
ylabel(ax2,'$\mathcal{E}_u-\mathcal{E}$','Interpreter','Latex','FontSize', 14)
xlabel(ax2,'$\mathcal{I}_u-\mathcal{E}$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(2,2,3);
%set(ax3,'YScale','log')
%set(ax3,'XScale','log')
hold(ax3,'on')
grid(ax3,'on')
ylabel(ax3,'$\Pi$','Interpreter','Latex','FontSize', 14)
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)

% ax4 = subplot(2,2,4);
% %set(ax4,'YScale','log')
% %set(ax4,'XScale','log')
% hold(ax4,'on')
% grid(ax4,'on')
% ylabel(ax4,'$\mathcal{I}_p-\mathcal{E}$','Interpreter','Latex','FontSize', 14)
% xlabel(ax4,'$t$','Interpreter','Latex','FontSize', 14)

fig2=figure;
set(fig2, 'Position',  [100,100,683,383])
set(gcf,'color','w');
axis1 = subplot(2,3,1);
%set(ax1,'YScale','log')
%set(ax1,'XScale','log')
hold(axis1,'on')
grid(axis1,'on')
ylabel(axis1,'$\mathcal{L}$','Interpreter','Latex','FontSize', 14)
xlabel(axis1,'$\dot{S}$','Interpreter','Latex','FontSize', 14)

axis2 = subplot(2,3,2);
%set(ax3,'YScale','log')
%set(ax3,'XScale','log')
hold(axis2,'on')
grid(axis2,'on')
ylabel(axis2,'$\mathcal{L}_u$','Interpreter','Latex','FontSize', 14)
xlabel(axis2,'$t$','Interpreter','Latex','FontSize', 14)

axis3 = subplot(2,3,3);
%set(ax3,'YScale','log')
%set(ax3,'XScale','log')
hold(axis3,'on')
grid(axis3,'on')
ylabel(axis3,'$\int_{0}^{50}\sqrt{\mathcal{I}_u(t)}dt$','Interpreter','Latex','FontSize', 14)
xlabel(axis3,'$t$','Interpreter','Latex','FontSize', 14)

axis4 = subplot(2,3,[4 5 6]);
%set(ax3,'YScale','log')
%set(ax3,'XScale','log')
hold(axis4,'on')
grid(axis4,'on')
ylabel(axis4,'$\Pi_u-\Pi$','Interpreter','Latex','FontSize', 14)
xlabel(axis4,'$t$','Interpreter','Latex','FontSize', 14)


a=0.1;
tx=6;
tspan = [0 50];   
%t=0:0.05:5;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
%tspan=logspace(0,5);
maxn=50;
%crls=rand([maxn,3]);

rng('default');
rng(1);
legaux={};
for n=1:maxn
    %A=[0,1;-1,-2];
    %A=[0,0.5,0.5; 0,0,0.5; 0,1,0];
    %positive definite random matrix
    % il=rand(n);
    % A=il*il.';
    %positive definite hurwitz matrix
    csys=rss(n);
    %csys = canon(sys,'companion');
    A=csys.A;
    B=csys.B; ninputs=size(B);
    %u=ones(ninputs(2),1);
    u=zeros(ninputs(2),1);

    Sigma0=0.1*eye(n);
    %X0=-5 +10*rand(1,n);
    X0=ones(1,n);
    At=A';
    D=0.01*eye(n);
    y0=Sigma0(:)';
    [t,y] = ode45(@(t,y) Sigma(t,y,A,D,n), tspan, y0, opts);
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

    for k=1:length(dSdt)
        Sx=reshape(y(k,:),n,n);
        dSdt(k)=0.5*trace(Sx\(A*Sx+Sx*A'+2*D));
        Pi(k)=((A*M(k,:)')'/D)*(A*M(k,:)')+trace(At/D*A*Sx+Sx\D+2*A);
        %Pit(k)=trace((A*M(k,:))*(A*M(k,:))')*trace(inv(D))+trace(inv(Sx))*trace(Sx*A*Sx*A'+2*Sx*A*D+D^2)*trace(inv(D));
        Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+(1/4)*trace(inv(Sx))*trace((A*Sx+Sx*At+2*D)^2)*trace(inv(D));
        Phi(k)=((A*M(k,:)')'/D)*(A*M(k,:)')+trace(At/D*A*Sx+A);
        E(k)=(A*M(k,:)')'/Sx*(A*M(k,:)')+0.5*trace((Sx\(A*Sx+Sx*A'+2*D))^2);
        Ip(k)=(1/n)*(trace(inv(Sx))*Pit(k)*trace(D))+dSdt(k)^2;
        Ep(k)=(1/n)*(trace(inv(Sx))*Pi(k)*trace(D))+dSdt(k)^2-2*S2(Sx\D+A);
    end
    
%     plot(ax1,dSdt,E,'Color',[crls(n,:)])    
%     plot(ax2,t,Ep-E,'Color',[crls(n,:)])    
%     plot(ax3,t,Pi,'Color',[crls(n,:)])
%     plot(ax4,t,Ip-E,'Color',[crls(n,:)])

    plot(ax1,dSdt,E,'k')
    plot(ax1,dSdt(end),E(end),'b*')
    plot(ax1,dSdt(1),E(1),'r*')
    plot(ax2,Ip-E,Ep-E,'k')   
    plot(ax2,Ip(end)-E(end),Ep(end)-E(end),'b*')
    plot(ax2,Ip(1)-E(1),Ep(1)-E(1),'r*')  
    plot(ax3,t,Pi,'k')
    
    plot(axis1,t,cumtrapz(sqrt(E)),'k')
    plot(axis2,t,cumtrapz(sqrt(Ep)),'k')
    plot(axis3,t,cumtrapz(sqrt(Ip)),'k')
    plot(axis4,t,Pit-Pi,'k')
    %plot(ax4,Ip-E,Ep-E,'k')
    
    %legaux{n-1}=sprintf('n=%d', n);
    
end

%sprintf('Some text %g', someValue)
    %leg1 = legend(legaux);
    %set(leg1,'Location','southwest');
    %leg1.Location = 'northeastoutside';
function dydt = Sigma(t,y,A,D,n)
   At=A';
   aux2=reshape(y,n,n);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*D(:);
end
 
function dydt = Mu(t,y,A,B,u,a,tx)
   dydt=A*y+B*(u*(1/(abs(a)*sqrt(pi))).*exp(-((t-tx)./a).^2));
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
% 
% function dydt = Mu(t,y,A,n)
%     dydt=zeros(n,length(t));
%     for k=1:length(t)
%      dydt(:,k)=expm(A*t(k))*y';
%     end
% end
