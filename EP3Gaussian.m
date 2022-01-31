%Entropy production three dimensional optical trap
%Author Adrian Guel November 2021

clear all;
close all;
clc
a=-10;
b=-3;
c=-1;
n=3;
Sigma0=diag([0.1 0.1 0.1]);
X0=[1,0.1,0.5];
A=diag([a b c]);
At=A';
D=diag([0.01 0.01 0.01]);
y0=Sigma0(:)';
tspan = [0 10];
[t,y] = ode45(@(t,y) Sigma(t,y,A,D,n), tspan, y0);
y(:,[2 3 4 6 7 8])=[];
[t,M] = ode45(@(t,y) Mu(t,y,A), t, X0);
fig=figure;
set(fig, 'Position',  [100,100,683,250])
set(gcf,'color','w');

dSdt=zeros(length(y(:,1)),1);
Pi=zeros(length(y(:,1)),1);
Phi=zeros(length(y(:,1)),1);
E=zeros(length(y(:,1)),1);
Pit=zeros(length(y(:,1)),1);
Ip=zeros(length(y(:,1)),1);
Ep=zeros(length(y(:,1)),1);

dSkdt=zeros(length(y(:,1)),n);
Pik=zeros(length(y(:,1)),n);
Eaux=zeros(length(y(:,1)),n);
Piksum=zeros(length(y(:,1)),1);
DSkdtsum=zeros(length(y(:,1)),1);
Esumi=zeros(length(y(:,1)),1);

for k=1:length(dSdt)
    Sx=diag([y(k,1) y(k,2) y(k,3)]);
    dSdt(k)=0.5*trace(Sx\(2*A*Sx+2*D));
    Pi(k)=((A*M(k,:)')'/D)*(A*M(k,:)')+trace(At/D*A*Sx+Sx\D+2*A);
    %Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+trace(inv(Sx))*trace(Sx*A*Sx*A'+2*Sx*A*D+D^2)*trace(inv(D));
    Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+(1/4)*trace(inv(Sx))*((trace(2*A*Sx+2*D))^2)*trace(inv(D));
    Phi(k)=((A*M(k,:)')'/D)*(A*M(k,:)')+trace(At/D*A*Sx+A);
    E(k)=(A*M(k,:)')'*(Sx\(A*M(k,:)'))+0.5*trace((Sx\(2*A*Sx+2*D))^2);
    Ip(k)=(1/n)*(trace(inv(Sx))*Pit(k)*trace(D))+dSdt(k)^2;
    Ep(k)=(1/n)*(trace(inv(Sx))*Pi(k)*trace(D))+dSdt(k)^2-2*S2(Sx\D+A);
end


for i=1:n
    dSkdt(:,i)=0.25*((2*A(i,i)*(y(:,i))+2*D(i,i))./(y(:,i))).^2;
    Pik(:,i)=((A(i,i)*M(:,i)).^2)/D(i,i)+(1/4)*((2*A(i,i)*(y(:,i))+2*D(i,i)).^2)./((y(:,i))*D(i,i));
    Eaux(:,i)=((A(i,i)*M(:,i)).^2)./(y(:,i))+0.5*((2*A(i,i)*(y(:,i))+2*D(i,i))./(y(:,i))).^2;
    Piksum=Piksum+(Pik(:,i)*D(i,i))./(y(:,i));
    DSkdtsum=DSkdtsum+dSkdt(:,i);
    Esumi=Esumi+Eaux(:,i);
end
Esumd=DSkdtsum+Piksum;

subplot(1,2,1)
colororder({'k','b'})
yyaxis left
plot(t,dSdt,'.')
hold on
plot(t,Pi-Phi,'k-')
ylabel('$\dot{S} (t)$','Interpreter','Latex','FontSize', 14)
yyaxis right
plot(t,E,'r.',t,Esumd,'b-')
ylabel('$\mathcal{E}(t)$','Interpreter','Latex','FontSize', 14)
xlabel('$t$','Interpreter','Latex','FontSize', 14)
leg1 = legend('$\dot{\mathcal{S}}(t)$','$\Pi-\Phi$','$\mathcal{E}$','$\bar{\mathcal{E}}$');
set(leg1,'Interpreter','latex');
grid on

subplot(1,2,2)
%yyaxis left
plot(t,Ep-Esumd)
%ylim([-1 1])
ylabel('$\mathcal{E}_u-\bar{\mathcal{E}}$','Interpreter','Latex','FontSize', 14)
%yyaxis right
%plot(t,Ip-Esumd)
%ylim([-1 1])
%ylabel('$\mathcal{I}_p-\mathcal{E}_m$','Interpreter','Latex','FontSize', 14)
xlabel('$t$','Interpreter','Latex','FontSize', 14)
%leg1 = legend('$\mathcal{E}_p-\mathcal{E}_m$','$\mathcal{I}_p-\mathcal{E}_m$');
%leg1 = legend('$\mathcal{E}_p-\mathcal{E}_m$');
%set(leg1,'Interpreter','latex');
grid on
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = (t < 8) & (t > 6); % range of t near perturbation
plot(t(indexOfInterest),Ep(indexOfInterest)-Esumd(indexOfInterest)) % plot on new axes
axis tight

function dydt = Sigma(t,y,A,D,n)
   At=A';
   aux2=reshape(y,3,3);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*D(:);
end
 
function dydt = Mu(t,y,A)
   dydt=A*y;
end
