clc
clear all

N=25;
sw=1;
a=1;
b=1;
A=[cos(pi/12) -sin(pi/12); sin(pi/12) cos(pi/12)];
C1=[1 0];
C2=[0 1];
R=sw^2;
Ry=a^2;
Rz=b^2;
P_y_n(:,:,1)=1e6*eye(2); 
P_z_n(:,:,1)=1e6*eye(2); 
x1(:,1)=[10*cos(pi/12);10*sin(pi/12)]; 
x2(:,1)=[10*cos(pi/12);10*sin(pi/12)]; 

for k=1:N   
    Y_original(k) = 10*cos(pi*(k-1)/12);
    Z_original(k) = 10*sin(pi*(k-1)/12);
end
X=[Y_original;Z_original];


y=[7.1165 9.6022 8.9144 9.2717 6.3400 4.0484 0.3411 -0.6784 -5.7726 ,...
    -5.4925 -9.4523 -9.7232 -9.5054 -9.7908 -7.7300 -5.9779 -4.5535 -1.5042 ,...
    -0.7044 3.2406 8.3029 6.1925 9.1178 9.0904 9.0662]; 
z=[0.000 3.1398 6.3739 9.5877 10.1450 10.1919 9.0683 10.2254 7.5799 ,...
    7.7231 5.4721 3.3990 0.9172 -1.3551 -5.2708 -9.7011 -9.4256 -9.3053 ,...
    -9.3815 -9.8822 -8.1876 -8.7501 -4.5653 -1.9179 -1.0000];


for k=1:N
    Kf_y(:,k)=P_y_n(:,:,k)*C1'*(C1*P_y_n(:,:,k)*C1'+Ry)^(-1); 
    Ky(:,k)=(A*P_y_n(:,:,k)*C1')*(C1*P_y_n(:,:,k)*C1'+Ry)^(-1); 
    Py(:,:,k)=P_y_n(:,:,k)-(P_y_n(:,:,k)*C1')*(C1*P_y_n(:,:,k)*C1'+Ry)^(-1)*C1*P_y_n(:,:,k);
    P_y_n(:,:,k+1)=A*P_y_n(:,:,k)*A'+R-Ky(:,k)*(C1*P_y_n(:,:,k)*C1'+Ry)*Ky(:,k)'; 
end

for k=1:N
    Kf_z(:,k)=P_z_n(:,:,k)*C2'*(C2*P_z_n(:,:,k)*C2'+Rz)^(-1); 
    Kz(:,k)=(A*P_z_n(:,:,k)*C2')*(C2*P_z_n(:,:,k)*C2'+Rz)^(-1); 
    Pz(:,:,k)=P_z_n(:,:,k)-(P_z_n(:,:,k)*C2')*(C2*P_z_n(:,:,k)*C2'+Rz)^(-1)*C2*P_z_n(:,:,k); 
    P_z_n(:,:,k+1)=A*P_z_n(:,:,k)*A'+R-Kz(:,k)*(C2*P_z_n(:,:,k)*C2'+Rz)*Kz(:,k)'; 
end

for k=1:N    
    x_y_n(:,k)=x1(:,k)+Kf_y(:,k)*(y(k)-C1*x1(:,k)); 
    x1(:,k+1)=A*x1(:,k)+Ky(:,k)*(y(k)-C1*x1(:,k)); 
    x_z_n(:,k)=x2(:,k)+Kf_z(:,k)*(z(k)-C2*x2(:,k)); 
    x2(:,k+1)=A*x2(:,k)+Kz(:,k)*(z(k)-C2*x2(:,k)); 
end

figure(1); 
plot(Y_original,Z_original,'bo'),hold on;
scatter(y,z,'go','filled'),hold on;
plot(x_y_n(1,:),x_z_n(2,:),'x-','color','r'),hold on;
xlabel('y')
ylabel('z')
axis([-15,15,-15,15]),axis('equal')
grid on
h=legend({'original positions','measured position','estimated positions'});
title('result');


figure(2); 
subplot(6,1,1)
plot([0:N-1],Kf_y(1,:)),hold on; 
xlabel('k')
t=title('$K_{fy1}(k)$');
set(t,'Interpreter','latex')
subplot(6,1,2)
plot([0:N-1],Kf_y(2,:)),hold on; 
xlabel('k')
t=title('$K_{fy2}(k)$');
set(t,'Interpreter','latex')


subplot(6,1,3)
plot([0:N-1],Kf_z(1,:)),hold on; %plot Kfz(k)
xlabel('k')
t=title('$K_{fz1}(k)$');
set(t,'Interpreter','latex')

subplot(6,1,4)
plot([0:N-1],Kf_z(2,:)),hold on; %plot Kfz(k)
xlabel('k')
t=title('$K_{fz2}(k)$');
set(t,'Interpreter','latex')


subplot(6,1,5)
plot([0:N-1],Ky(1,:)),hold on; 
xlabel('k')
t=title('$K_{y1}(k)$');
set(t,'Interpreter','latex')


subplot(6,1,6)
plot([0:N-1],Ky(2,:)),hold on; 
xlabel('k')
t=title('$K_{y2}(k)$');
set(t,'Interpreter','latex')



figure(4); 
subplot(4,1,1)
plot([0:N-1],Kz(1,:)),hold on; 
xlabel('k')
t=title('$K_{z1}(k)$');
set(t,'Interpreter','latex')

subplot(4,1,2)
plot([0:N-1],Kz(2,:)),hold on; 
xlabel('k')
t=title('$K_{z2}(k)$');
set(t,'Interpreter','latex')

subplot(4,1,3)
Py1=Py(1,1,:);
Py2=Py(2,2,:);
plot([0:N-1],Py1(:)),hold on;
xlabel('k')
t=title('$P_{y11}(k|k)$');
set(t,'Interpreter','latex')

subplot(4,1,4)
plot([0:N-1],Py2(:)),hold on; 
xlabel('k')
t=title('$P_{y22}(k|k)$');
set(t,'Interpreter','latex')



figure(5);
subplot(6,1,1)
Pz1=Pz(1,1,:);
Pz2=Pz(2,2,:);
plot([0:N-1],Pz1(:)),hold on; 
xlabel('k')
t=title('$P_{z11}(k|k)$');
set(t,'Interpreter','latex')

subplot(6,1,2)
plot([0:N-1],Pz2(:)),hold on; 
xlabel('k')
t=title('$P_{z22}(k|k)$');
set(t,'Interpreter','latex')


subplot(6,1,3)
Pmy1=P_y_n(1,1,:);
Pmy2=P_y_n(2,2,:);
plot([0:N],Pmy1(:)),hold on; 
xlabel('k')
t=title('$P_{y11}(k+1|k)$');
set(t,'Interpreter','latex')


subplot(6,1,4)
plot([0:N],Pmy2(:)),hold on; 
xlabel('k')
t=title('$P_{y22}(k+1|k)$');
set(t,'Interpreter','latex')



subplot(6,1,5)
Pmz1=P_z_n(1,1,:);
Pmz2=P_z_n(2,2,:);
plot([0:N],Pmz1(:)),hold on; 
xlabel('k')
t=title('$P_{z11}(k+1|k)$');
set(t,'Interpreter','latex')


subplot(6,1,6)
plot([0:N],Pmz2(:)),hold on; 
xlabel('k')
t=title('$P_{z22}(k+1|k)$');
set(t,'Interpreter','latex');

Bias_y= sum(X(1,:)-x_y_n(1,:))./N
Bias_z= sum(X(2,:)-x_z_n(2,:))./N
Var_y = sum((X(1,:)-x_y_n(1,:)).^2)./N
Var_z = sum((X(2,:)-x_z_n(2,:)).^2)./N
