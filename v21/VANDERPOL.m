
function VANDERPOL

%
% Выбор метода решения:
method_type='SDIRC'; % Допустимые значения:
                        % 1) 'MidPoint' --- неявный метод средней точки
                        %      (одношаговый, одностадийный, 2-ого порядка)
                        % 2) 'HH' --- метод Хаммера-Холлингсворта
                        %      (одношаговый, двухстадийный, 4-ого порядка)
                        % 3) 'KB' --- метод Кунцмана-Бутчера
                        %      (одношаговый, трёхстадийный, 6-ого порядка)
                        % 4) 'SDIRC' --- диагональный неявный метод Рунге-Кутты
                        %      (одношаговый, двухстадийный, 2-ого порядка,
                        %       ассимптотически устойчивый)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Параметры уравнения:
x_init=2; % начальные условия
y_init=0;
a=1e3; % "параметр жёсткости"
T=20; % конец временного интервала
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Параметры численного решения:
N=5000;     % кол-во узлов равномерной сетки на отрезке [0,T]
tau=T/N;    % шаг сетки
toler=1e-9; % заданная относительная точность решения нелинейных уравнений 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Решение:
t=linspace(0,T,N+1);
Sol=[[x_init;y_init],zeros(2,N)];
%
if(strcmp(method_type,'MidPoint'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK1(Sol(1,i),Sol(2,i),a,tau,k,...
                                                 0.5),...
                             [0;0],optimset('Display','off','TolX',toler));
    
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+tau*k;
    end
    toc
end
%
if(strcmp(method_type,'HH'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK2(Sol(1,i),Sol(2,i),a,tau,k,...
                                    1/4,      1/4-sqrt(3)/6,...
                                1/4+sqrt(3)/6,      1/4   ), ...
                    [0;0;0;0],optimset('Display','off','TolX',toler));
    
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+tau*k(1:2)/2+tau*k(3:4)/2;
    end
    toc
end
%
if(strcmp(method_type,'SDIRC'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK2(Sol(1,i),Sol(2,i),a,tau,k,...
                                  1+sqrt(2)/2,            0      ,...
                                   -sqrt(2),      1+sqrt(2)/2   ),...
                             [0;0;0;0],optimset('Display','off','TolX',toler));
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+tau*k(1:2)/2+tau*k(3:4)/2;
    
    end
    toc
end
%
if(strcmp(method_type,'KB'))
    tic
    for i=1:N
    [k,~,exitflag,~]=fsolve(@(k) RK3(Sol(1,i),Sol(2,i),a,tau,k,...
                         5/36,        2/9-sqrt(15)/15, 5/36-sqrt(15)/30,...
                    5/36+sqrt(15)/24,      2/9,        5/36-sqrt(15)/24,...
                    5/36+sqrt(15)/30, 2/9+sqrt(15)/15,      5/36      ),...
                    [0;0;0;0;0;0],optimset('Display','off','TolX',toler));
    
    if(exitflag~=1)
        exitflag
    end
    Sol(:,i+1)=Sol(:,i)+5*tau*k(1:2)/18+4*tau*k(3:4)/9+5*tau*k(5:6)/18;
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
title('$x(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
%
figure(2)
hold on
plot(t,Sol(2,:),'b-')
title('$y(t)$','Interpreter','latex')
xlabel('$t$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
%
figure(3)
hold on
plot(Sol(1,:),Sol(2,:),'b.')
plot(x_init,y_init,'r+')
title('Trajectory')
xlabel('$x(t)$','Interpreter','latex')
ylabel('$y(t)$','Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [F]=RK1(x,y,a,tau,k,c1)
%
x=x+tau*c1*k(1);
y=y+tau*c1*k(2);
F=[-a*(x^3/3-x)+y-k(1);
   -a*x-k(2)];
%
function [F]=RK2(x,y,a,tau,k,c11,c12,c21,c22)
%
x1=x+c11*tau*k(1)+c12*tau*k(3);
y1=y+c11*tau*k(2)+c12*tau*k(4);
x2=x+c21*tau*k(1)+c22*tau*k(3);
y2=y+c21*tau*k(2)+c22*tau*k(4);
F=[-a*(x1^3/3-x1)+y1-k(1);
   -a*x1-k(2);
   -a*(x2^3/3-x2)+y2-k(3);
   -a*x2-k(4)];
%
function [F]=RK3(x,y,a,tau,k,c11,c12,c13,c21,c22,c23,c31,c32,c33)
%
x1=x+c11*tau*k(1)+c12*tau*k(3)+c13*tau*k(5);
y1=y+c11*tau*k(2)+c12*tau*k(4)+c13*tau*k(6);
x2=x+c21*tau*k(1)+c22*tau*k(3)+c23*tau*k(5);
y2=y+c21*tau*k(2)+c22*tau*k(4)+c23*tau*k(6);
x3=x+c31*tau*k(1)+c32*tau*k(3)+c33*tau*k(5);
y3=y+c31*tau*k(2)+c32*tau*k(4)+c33*tau*k(6);
F=[-a*(x1^3/3-x1)+y1-k(1);
   -a*x1-k(2);
   -a*(x2^3/3-x2)+y2-k(3);
   -a*x2-k(4);
   -a*(x3^3/3-x3)+y3-k(5);
   -a*x3-k(6)];