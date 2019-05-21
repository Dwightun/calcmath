function VANDERPOL
%
% Выбор метода решения:
method_type='Rado5'; % Допустимые значения:
                        % 1) 'Euler' --- неявный метод Эйлер
                        %      (одношаговый, одностадийный, 2-ого порядка)
                        % 2) 'Rado3' --- метод Радо IIA
                        %      (одношаговый, двухстадийный, 3-ого порядка)
                        % 3) 'Rado5' --- метод Радо IIA
                        %      (одношаговый, трёхстадийный, 5-ого порядка)
                        % 4) 'SDIRC' --- однакратно диагональный неявный метод Рунге-Кутты
                        %      (одношаговый, двухстадийный, 2-ого порядка)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Параметры уравнения:
x_init=1; % начальные условия
y_init=0;
z_init=1;
w_init=0;
a=102;      % "параметр жёсткости"
b=100;      % "параметр жёсткости"
k12=0.5;
k21=0.5;
k22=1.0;
T=180; % конец временного интервала
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Параметры численного решения:
N=5000;     % кол-во узлов равномерной сетки на отрезке [0,T]
tau=T/N;    % шаг сетки
toler=1e-5; % заданная относительная точность решения нелинейных уравнений 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Решение:
t=linspace(0,T,N+1);
Sol=[[x_init;y_init;z_init;w_init],zeros(4,N)];
%
if(strcmp(method_type,'Euler'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK1(Sol(1,i),Sol(2,i),Sol(3,i),Sol(4,i),...
                                            a,b,k12,k21,k22,tau,k,...
                                                 1),...
                             zeros(4,1),optimset('Display','off','TolX',toler));
    i
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+tau*k;
    end
    toc
end
%
if(strcmp(method_type,'Rado3'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK2(Sol(1,i),Sol(2,i),Sol(3,i),Sol(4,i),...
                                        a,b,k12,k21,k22,tau,k,...
                                          5/12,      -1/12,...
                                          3/4,        1/4   ), ...
                    zeros(8,1),optimset('Display','off','TolX',toler));
    i
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+3*tau*k(1:4)/4+tau*k(5:8)/4;
    end
    toc
end
%
if(strcmp(method_type,'SDIRC'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK2(Sol(1,i),Sol(2,i),Sol(3,i),Sol(4,i),...
                                            a,b,k12,k21,k22,tau,k,...
                                     1-sqrt(2)/2,            0      ,...
                                        sqrt(2),      1-sqrt(2)/2   ),...
                             zeros(8,1),optimset('Display','off','TolX',toler));
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+tau*k(1:4)/2+tau*k(5:8)/2;
    i
    end
    toc
end
%
if(strcmp(method_type,'Rado5'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK3(Sol(1,i),Sol(2,i),Sol(3,i),Sol(4,i),...
                                        a,b,k12,k21,k22,tau,k,...
     (88-7*sqrt(6))/360,      (296-169*sqrt(6))/1800,    (-2+3*sqrt(6))/225,...
   (296+169*sqrt(6))/1800,      (88+7*sqrt(6))/360,      (-2-3*sqrt(6))/225,...
      (16-sqrt(6))/36,            (16+sqrt(6))/36,               1/9      ),...
                    zeros(12,1),optimset('Display','off','TolX',toler));
    i
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+(16-sqrt(6))*tau*k(1:4)/36+(16+sqrt(6))*tau*k(5:8)/36+tau*k(9:12)/9;
    end
    toc
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Отрисовка результатов:
figure(1)
hold on
plot(t,Sol(1,:),'b-')
plot(t,Sol(3,:),'r-')
title('$x(t),z(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$x,z$','Interpreter','latex')
%
figure(2)
hold on
plot(t,Sol(2,:),'b-')
plot(t,Sol(4,:),'r-')
title('$y(t),w(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$y,w$','Interpreter','latex')
%
figure(3)
hold on
plot(Sol(1,:),Sol(2,:),'b.')
plot(x_init,y_init,'k+')
title('Траектория в плоскости (x,y)')
xlabel('$x(t)$','Interpreter','latex')
ylabel('$y(t)$','Interpreter','latex')
%
figure(4)
hold on
plot(Sol(3,:),Sol(4,:),'r.')
plot(x_init,y_init,'k+')
title('Траектория в плоскости (z,w)')
xlabel('$z(t)$','Interpreter','latex')
ylabel('$w(t)$','Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [F]=RK1(x,y,z,w,a,b,k12,k21,k22,tau,k,c1)
%
x=x+tau*c1*k(1);
y=y+tau*c1*k(2);
z=z+tau*c1*k(3);
w=w+tau*c1*k(4);
F=[-a*(x^3/3-x)+y-k(1);
      -x-k12*z-k(2);
   -b*(z^3/3-z)+w-k(3);
     -k22*z-k21*x-k(4);];
%
function [F]=RK2(x,y,z,w,a,b,k12,k21,k22,tau,k,c11,c12,c21,c22)
%
x1=x+c11*tau*k(1)+c12*tau*k(5);
y1=y+c11*tau*k(2)+c12*tau*k(6);
z1=z+c11*tau*k(3)+c12*tau*k(7);
w1=w+c11*tau*k(4)+c12*tau*k(8);
x2=x+c21*tau*k(1)+c22*tau*k(5);
y2=y+c21*tau*k(2)+c22*tau*k(6);
z2=z+c21*tau*k(3)+c22*tau*k(7);
w2=w+c21*tau*k(4)+c22*tau*k(8);
F=[-a*(x1^3/3-x1)+y1-k(1);
   -x1-k12*z1-k(2);
   -b*(z1^3/3-z1)+w1-k(3);
   -k22*z1-k21*x1-k(4);
   -a*(x2^3/3-x2)+y2-k(5);
   -x2-k12*z2-k(6);
   -b*(z2^3/3-z2)+w2-k(7);
   -k22*z2-k21*x2-k(8);];
%
function [F]=RK3(x,y,z,w,a,b,k12,k21,k22,tau,k,c11,c12,c13,c21,c22,c23,c31,c32,c33)
%
x1=x+c11*tau*k(1)+c12*tau*k(5)+c13*tau*k(9);
y1=y+c11*tau*k(2)+c12*tau*k(6)+c13*tau*k(10);
z1=z+c11*tau*k(3)+c12*tau*k(7)+c13*tau*k(11);
w1=w+c11*tau*k(4)+c12*tau*k(8)+c13*tau*k(12);
x2=x+c21*tau*k(1)+c22*tau*k(5)+c23*tau*k(9);
y2=y+c21*tau*k(2)+c22*tau*k(6)+c23*tau*k(10);
z2=z+c21*tau*k(3)+c22*tau*k(7)+c23*tau*k(11);
w2=w+c21*tau*k(4)+c22*tau*k(8)+c23*tau*k(12);
x3=x+c31*tau*k(1)+c32*tau*k(5)+c33*tau*k(9);
y3=y+c31*tau*k(2)+c32*tau*k(6)+c33*tau*k(10);
z3=z+c31*tau*k(3)+c32*tau*k(7)+c33*tau*k(11);
w3=w+c31*tau*k(4)+c32*tau*k(8)+c33*tau*k(12);
F=[-a*(x1^3/3-x1)+y1-k(1);
   -x1-k12*z1-k(2);
   -b*(z1^3/3-z1)+w1-k(3);
   -k22*z1-k21*x1-k(4);
   -a*(x2^3/3-x2)+y2-k(5);
   -x2-k12*z2-k(6);
   -b*(z2^3/3-z2)+w2-k(7);
   -k22*z2-k21*x2-k(8);
    -a*(x3^3/3-x3)+y3-k(9);
   -x3-k12*z3-k(10);
   -b*(z3^3/3-z3)+w3-k(11);
   -k22*z3-k21*x3-k(12);];