clc
clear all

N=10000; 
n=1; 
T=1; 
sw=0.1; 
sv=1; 
A=[1 T;0 1]; 
C=[1 0];
R1=sw^2*[T^4/4 T^3/2;T^3/2 T^2]; 
R2=sv^2; 
P_n(:,:,1)=1e6*eye(2); 
w=random('normal',0,sw,n,N); 
v=random('normal',0,sv,n,N); 

for k=1:N
    Kf(:,k)=P_n(:,:,k)*C'*(C*P_n(:,:,k)*C'+R2)^(-1); 
    K(:,k)=(A*P_n(:,:,k)*C')*(C*P_n(:,:,k)*C'+R2)^(-1); 
    P(:,:,k)=P_n(:,:,k)-(P_n(:,:,k)*C')*(C*P_n(:,:,k)*C'+R2)^(-1)*C*P_n(:,:,k); 
    P_n(:,:,k+1)=A*P_n(:,:,k)*A'+R1-K(:,k)*(C*P_n(:,:,k)*C'+R2)*K(:,k)'; 
end

for run=1:n
    x_m(:,1)=[0 0]'; 
    x(:,1)=[0 30]'; 
    for k=1:N
        x(:,k+1)=A*x(:,k)+[T^2/2 T]'*w(run,k); 
        y(k)=C*x(:,k)+v(run,k); 
        xh(:,k)=x_m(:,k)+Kf(:,k)*(y(k)-C*x_m(:,k)); 
        x_m(:,k+1)=A*x_m(:,k)+K(:,k)*(y(k)-C*x_m(:,k)); 
    end
end


figure();
subplot(2, 1, 1);
plot([0:N-1],x(1,1:N),'b-'),hold on; 
plot([0:N-1],xh(1,1:N),'r:'),hold on;
xlabel('k')
h=legend({'$x_{1}(k)$','$\hat{x}_{1}(k|k)$'});
set(h,'Interpreter','latex')
t=title('Graph 4');
set(t,'Interpreter','latex')

subplot(2, 1, 2);
plot([0:N-1],x(2,1:N),'b-'),hold on; 
plot([0:N-1],xh(2,1:N),'r:'),hold on;
xlabel('k')
h=legend({'$x_{2}(k)$','$\hat{x}_{2}(k|k)$'});
set(h,'Interpreter','latex')
t=title('Graph 5');
set(t,'Interpreter','latex')

figure();
subplot(2, 1, 1);
P1=P(1,1,:);
P2=P(2,2,:);
plot([0:N-1],P1(:),'bo'),hold on; 
xlabel('k')
h=legend({'$P_{11}(k|k)$'});
t=title('Graph 6');
set(t,'Interpreter','latex')
set(h,'Interpreter','latex')

subplot(2, 1, 2);
plot([0:N-1],P2(:),'bo'),hold on; 
xlabel('k')
h=legend({'$P_{22}(k|k)$'});
t=title('Graph 7');
set(t,'Interpreter','latex')
set(h,'Interpreter','latex')

figure();
subplot(2, 1, 1);
plot([0:N-1],Kf(1,:),'bo'),hold on; 
xlabel('k')
h=legend({'$K_{f1}(k)$'});
t=title('Graph 8');
set(t,'Interpreter','latex')
set(h,'Interpreter','latex')

subplot(2, 1, 2);
plot([0:N-1],Kf(2,:),'bo')
hold on; %plot Kf(k)
xlabel('k')
h=legend({'$K_{f2}(k)$'});
t=title('Graph 9');
set(t,'Interpreter','latex')
set(h,'Interpreter','latex')

Bias_1= sum(x(1,1:N)-xh(1,1:N))./N
Bias_2= sum(x(2,1:N)-xh(2,1:N))./N
Var_1 = sum((x(1,1:N)-xh(1,1:N)).^2)./N
Var_2 = sum((x(2,1:N)-xh(2,1:N)).^2)./N
