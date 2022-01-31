%Entropy production brownian oscillator abrupt events
%Author Adrian Guel December 2021

clear all;
close all;
clc
a11=0;
a12=1;
a21=-1;
a22=-1;
n=2;
Sigma0=[0.1,0;0,0.1];
X0=[1,1];
A=[a11,a12;a21,a22];
B=[0;0];
At=A';
D=[0.01,0;0,0.01];
a=0.1;
b=[0,0];
tx=[1,1];
y0=Sigma0(:)';
tspan = [0 8];
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
[t,y] = ode45(@(t,y) Sigma(t,y,A,D,n,a,b,tx), tspan, y0,opts);

fig=figure;
set(fig, 'Position',  [100,100,700,283])
set(gcf,'color','w');
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
% 
% subplot(2,2,2)
 [t,M] = ode45(@(t,y) Mu(t,y,A,B,a,tx), t, X0);
% yyaxis left
% plot(t,M(:,1),'--')
% ylabel('$\mu_1(t)$','Interpreter','Latex','FontSize', 14)
% yyaxis right
% plot(t,M(:,2),'-.')
% ylabel('$\mu_2(t)$','Interpreter','Latex','FontSize', 14)
% xlabel('$t$','Interpreter','Latex','FontSize', 14)
% leg1 = legend('$\mu_1$','$\mu_{2}$');
% set(leg1,'Interpreter','latex');
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

Dx1=D(1,1)+(b(1)/(abs(a)*sqrt(pi))).*exp(-((t-tx(1))./a).^2);
Dx2=D(2,2)+(b(2)/(abs(a)*sqrt(pi))).*exp(-((t-tx(2))./a).^2);

subplot(1,2,1)
yyaxis left
plot(t,dSdt,'-.',t,Dx2,'k-')
ylabel('$\dot{S} (t),D_{22}(t)$','Interpreter','Latex','FontSize', 14)
yyaxis right
plot(t,Pi,'-.',t,Phi,'b--')
ylabel('$\Pi (t),\Phi(t)$','Interpreter','Latex','FontSize', 14)
xlabel('$t$','Interpreter','Latex','FontSize', 14)
leg1 = legend('$\dot{S}$','$D_{22}$','$\Pi$','$\Phi$');
set(leg1,'Interpreter','latex');
grid on

subplot(1,2,2)
yyaxis left
plot(t,Dx2,'k-')
ylabel('$D_{22}(t)$','Interpreter','Latex','FontSize', 14)
yyaxis right
plot(t,E,'b-',t,Sp,'.')
ylabel('$\mathcal{E} (t)$','Interpreter','Latex','FontSize', 14)
xlabel('$t$','Interpreter','Latex','FontSize', 14)
leg1 = legend('$D_{22} (t)$','$\mathcal{E}(t)$','$\mathcal{S}_p(t)$');
set(leg1,'Interpreter','latex');
grid on

function dydt = Sigma(t,y,A,D,n,a,b,tx)
   At=A';
   aux2=reshape(y,2,2);
   Dx=Dfunction2(D,t,b,a,tx);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*Dx(:);
end
 
function dydt = Mu(t,y,A,B,a,tx)
   dydt=A*y+B*(1/(abs(a)*sqrt(pi))).*exp(-((t-tx(1))./a).^2);
end
