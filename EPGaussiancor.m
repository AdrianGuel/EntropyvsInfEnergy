%Entropy production brownian oscillator
%correlation analysis
%Author Adrian Guel November 2021

clear all;
close all;

fig=figure;
set(fig, 'Position',  [615,328,691,354])
set(gcf,'color','w');

ax1 = subplot(3,2,[2 4 6]);
%set(ax1,'YScale','log')
%set(ax1,'XScale','log')
hold(ax1,'on')
grid(ax1,'on')
xlabel(ax1,'$x_1\sim \mathcal{N}(\mu,\Sigma)$','Interpreter','Latex','FontSize', 14)
ylabel(ax1,'$x_2\sim \mathcal{N}(\mu,\Sigma)$','Interpreter','Latex','FontSize', 14)
axis(ax1,'square')

ax2 = subplot(3,2,1);
%set(ax2,'YScale','log')
%set(ax1,'XScale','log')
hold(ax2,'on')
grid(ax2,'on')
xlabel(ax2,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax2,'$|\mathcal{E}-\bar{\mathcal{E}}|$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(3,2,3);
%set(ax3,'YScale','log')
%set(ax1,'XScale','log')
hold(ax3,'on')
grid(ax3,'on')
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax3,'$\mathcal{C}_\Pi$','Interpreter','Latex','FontSize', 14)

ax4 = subplot(3,2,5);
%set(ax3,'YScale','log')
%set(ax1,'XScale','log')
hold(ax4,'on')
grid(ax4,'on')
xlabel(ax4,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax4,'$\mathcal{C}_S$','Interpreter','Latex','FontSize', 14)

% ax5 = subplot(3,3,[3 6 9]);
% %set(ax1,'YScale','log')
% %set(ax1,'XScale','log')
% hold(ax5,'on')
% grid(ax5,'on')
% xlabel(ax5,'$x_1$','Interpreter','Latex','FontSize', 14)
% ylabel(ax5,'$x_2$','Interpreter','Latex','FontSize', 14)
% zlabel(ax5,'$p(\mathbf{x};t)$','Interpreter','Latex','FontSize', 14)
% axis(ax5,'square')
% v = [-5 -2 5];
% view(v)

a11=0;
a12=1;
%a21=[0 -1 -2 -3];
a21=-2;
%a22=[0 -1 -2 -3];
a22=-1;
n=2;
Sigma0=[0.1,0;0,0.1];
X0=[1,1];
%A=[a11,a12;a21,a22];
D=[0.01,0;0,0.01];
y0=Sigma0(:)';


for l=1:length(a22)
    %A=[a11,a12;a21(l),a22];
    A=[a11,a12;a21,a22(l)];
    At=A';
    tspan = [0 50];
    [t,y] = ode45(@(t,y) Sigma(t,y,A,D,n), tspan, y0);
    [t,M] = ode45(@(t,y) Mu(t,y,A), t, X0);

    x1 = -2:0.1:2;
    x2 = -2:0.1:2;
    [X1,X2] = meshgrid(x1,x2);
    X = [X1(:) X2(:)];  
    for k=1:10:length(t)
        yx = mvnpdf(X,[M(k,1),M(k,2)],[y(k,1),y(k,2);y(k,3),y(k,4)]);
        yx = reshape(yx,length(x2),length(x1));
        contour(ax1,real(x1),real(x2),real(yx),10,'k')
        txt = ['t=' num2str(t(k))];
        text(ax1,M(k,1),M(k,2),txt, 'Color', 'r')
        %surf(ax5,real(x1),real(x2),real(yx))
        %hold on              
    end          
    plot(ax1,M(:,1),M(:,2),'b','LineWidth',2)

    dSdt=zeros(length(y(:,1)),1);
    Pi=zeros(length(y(:,1)),1);
    Phi=zeros(length(y(:,1)),1);
    E=zeros(length(y(:,1)),1);
    Pit=zeros(length(y(:,1)),1);
    Ep=zeros(length(y(:,1)),1);
    Sp=zeros(length(y(:,1)),1);
    
    bE=zeros(length(y(:,1)),1);
    bS=zeros(length(y(:,1)),1);
    bPi=zeros(length(y(:,1)),1);
    
    for k=1:length(t)
        Sx=[y(k,1),y(k,3);y(k,2),y(k,4)];
        dSdt(k)=0.5*trace(Sx\(A*Sx+Sx*A'+2*D));
        Pi(k)=trace((A*M(k,:)')'*(D\(A*M(k,:)'))+A'*(D\A)*Sx+Sx\D+2*A);
        %Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+trace(inv(Sx))*trace(Sx*A*Sx*A'+2*Sx*A*D+D^2)*trace(inv(D));
        Pit(k)=trace((A*M(k,:)')*(A*M(k,:)')')*trace(inv(D))+(1/4)*trace(inv(Sx))*((trace(A*Sx+Sx*A'+2*D))^2)*trace(inv(D));
        Phi(k)=trace((A*M(k,:)')'*(D\(A*M(k,:)'))+A'*(D\A)*Sx+A);
        E(k)=(A*M(k,:)')'*(Sx\(A*M(k,:)'))+0.5*trace((Sx\(A*Sx+Sx*A'+2*D))^2);
        Ep(k)=(1/n)*(trace(inv(Sx))*Pit(k)*trace(D))+dSdt(k)^2;
        Sp(k)=(1/n)*(trace(inv(Sx))*Pi(k)*trace(D))+dSdt(k)^2-2*S2(Sx\D+A);
       
            bE(k)=((A(1,1)*M(k,1)+A(1,2)*M(k,2))^2)/Sx(1,1)+0.5*((2*D(1,1)+2*A(1,1)*Sx(1,1)+A(1,2)*(Sx(1,2)+Sx(2,1)))/Sx(1,1))^2....
            +((A(2,1)*M(k,1)+A(2,2)*M(k,2))^2)/Sx(2,2)+0.5*((2*D(2,2)+2*A(2,2)*Sx(2,2)+A(2,1)*(Sx(1,2)+Sx(2,1)))/Sx(2,2))^2;
            bPi(k)=((A(1,1)*M(k,1)+A(1,2)*M(k,2))^2)/D(1,1)+0.25*((2*D(1,1)+2*A(1,1)*Sx(1,1)+A(1,2)*(Sx(1,2)+Sx(2,1))))^2/(Sx(1,1)*D(1,1))....
                +((A(2,1)*M(k,1)+A(2,2)*M(k,2))^2)/D(2,2)+0.25*((2*D(2,2)+2*A(2,2)*Sx(2,2)+A(2,1)*(Sx(1,2)+Sx(2,1))))^2/(Sx(2,2)*D(2,2));
            bS(k)=0.5*((2*D(1,1)+2*A(1,1)*Sx(1,1)+A(1,2)*(Sx(1,2)+Sx(2,1)))/Sx(1,1))....
                +0.5*((2*D(2,2)+2*A(2,2)*Sx(2,2)+A(2,1)*(Sx(1,2)+Sx(2,1)))/Sx(2,2));
            
    end
     plot(ax2,t,abs(E-bE),'b','LineWidth',2)
     plot(ax3,t,abs(Pi-bPi),'b','LineWidth',2)
     plot(ax4,t,abs(dSdt-bS),'b','LineWidth',2)

end


function dydt = Sigma(t,y,A,D,n)
   At=A';
   aux2=reshape(y,2,2);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*D(:);
end
 
function dydt = Mu(t,y,A)
   dydt=A*y;
end
