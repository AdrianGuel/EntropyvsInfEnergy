%Entropy production brownian oscillator
%Author Adrian Guel November 2021

clear all;
close all;
clc
a11=-1;
a12=1;
a21=-1;
a22=0;
%a22=[0 -1 -2 -3];
n=2;
Sigma0=[0.1,0;0,0.1];
X0=[1,1];
%A=[a11,a12;a21,a22];
D=[0.01,0;0,0.01];
y0=Sigma0(:)';


fig=figure;
set(fig, 'Position',  [100,100,683,383])
set(gcf,'color','w');

for l=1:length(a21)
    A=[a11,a12;a21(l),a22];
    At=A';
    tspan = [0 5];
    [t,y] = ode45(@(t,y) Sigma(t,y,A,D,n), tspan, y0);
    [t,M] = ode45(@(t,y) Mu(t,y,A), t, X0);


    % subplot(2,2,1)
    % yyaxis left
    % plot(t,y(:,1),'r-',t,y(:,4),'--')
    % ylabel('$\Sigma_{11}(t)$,$\Sigma_{22}(t)$','Interpreter','Latex','FontSize', 14)
    % yyaxis right
    % plot(t,y(:,2),'-.',t,y(:,3),'r-')
    % ylabel('$\Sigma_{12} (t)$','Interpreter','Latex','FontSize', 14)
    % xlabel('$t$','Interpreter','Latex','FontSize', 14)
    % leg1 = legend('$\Sigma_{11}$','$\Sigma_{12}$','$\Sigma_{21}$','$\Sigma_{22}$');
    % set(leg1,'Interpreter','latex');
    % grid on

    % subplot(1,2,2)
    %     x1 = -2:0.1:2;
    %     x2 = -2:0.1:2;
    %     [X1,X2] = meshgrid(x1,x2);
    %     X = [X1(:) X2(:)];  
    %     for k=1:1:length(t)
    %         yx = mvnpdf(X,[M(k,1),M(k,2)],[y(k,1),y(k,2);y(k,3),y(k,4)]);
    %         yx = reshape(yx,length(x2),length(x1));
    %         contour(real(x1),real(x2),real(yx),1,'k')
    %         hold on              
    %     end          
    % plot(M(:,1),M(:,2),'b','LineWidth',2)
    % % %xlabel('$t$','Interpreter','Latex','FontSize', 16)
    % xlabel('$x_1\sim \mathcal{N}(\mu,\Sigma)$','Interpreter','Latex','FontSize', 12)
    % ylabel('$x_2\sim \mathcal{N}(\mu,\Sigma)$','Interpreter','Latex','FontSize', 12)
    % %view(46.223544303797468,19.058734177215193)
    % %Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam) % a drawing function
    % grid on

    dSdt=zeros(length(y(:,1)),1);
    Pi=zeros(length(y(:,1)),1);
    Phi=zeros(length(y(:,1)),1);
    E=zeros(length(y(:,1)),1);
    Pit=zeros(length(y(:,1)),1);
    Ep=zeros(length(y(:,1)),1);
    Sp=zeros(length(y(:,1)),1);
    for k=1:length(dSdt)
        Sx=[y(k,1),y(k,3);y(k,2),y(k,4)];
        dSdt(k)=0.5*trace(Sx\(A*Sx+Sx*A'+2*D));
        Pi(k)=trace((A*M(k,:)')'*(D\(A*M(k,:)'))+A'*(D\A)*Sx+Sx\D+2*A);
        %Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+trace(inv(Sx))*trace(Sx*A*Sx*A'+2*Sx*A*D+D^2)*trace(inv(D));
        Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+(1/4)*trace(inv(Sx))*((trace(A*Sx+Sx*A'+2*D))^2)*trace(inv(D));
        Phi(k)=trace((A*M(k,:)')'*(D\(A*M(k,:)'))+A'*(D\A)*Sx+A);
        E(k)=(A*M(k,:)')'*(Sx\(A*M(k,:)'))+0.5*trace((Sx\(A*Sx+Sx*A'+2*D))^2);
        Ep(k)=(1/n)*(trace(inv(Sx))*Pit(k)*trace(D))+dSdt(k)^2;
        Sp(k)=(1/n)*(trace(inv(Sx))*Pi(k)*trace(D))+dSdt(k)^2-2*S2(Sx\D+A);
    end

    subplot(2,2,l)
    yyaxis left
    plot(t,dSdt,'k.')
    ylabel('$\dot{S} (t)$','Interpreter','Latex','FontSize', 14)
    hold on
    plot(t,Pi-Phi,'r-')
    %ylabel('$\Pi (t),\Phi(t)$','Interpreter','Latex','FontSize', 14)
    %xlabel('$t$','Interpreter','Latex','FontSize', 14)
    %leg1 = legend('$\dot{S}$','$\Pi-\Phi$');
    %set(leg1,'Interpreter','latex');
    %grid on

    %subplot(2,2,3)
    % yyaxis left
    % plot(t,Pi,'.')
    % ylabel('$\Pi (t)$','Interpreter','Latex','FontSize', 14)
     yyaxis right
    plot(t,Sp-E,'b-')
    text(5, 30, sprintf('a_{21}=%.2f', a21(l)))
    ylabel('$\mathcal{S}_p(t)-\mathcal{E} (t)$','Interpreter','Latex','FontSize', 14)
    grid on
end
    xlabel('$t$','Interpreter','Latex','FontSize', 14)
    leg1 = legend('$\dot{S}$','$\Pi-\Phi$','$\mathcal{S}_p(t)-\mathcal{E}(t)$');
    set(leg1,'Interpreter','latex');

function dydt = Sigma(t,y,A,D,n)
   At=A';
   aux2=reshape(y,2,2);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*D(:);
end
 
function dydt = Mu(t,y,A)
   dydt=A*y;
end
