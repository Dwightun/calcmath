function PROBLEM25
%
N=100;     % Число узлов сетки по x
h=1/N;    % Шаг сетки по x
x=linspace(0,1,N+1); % Узлы сетки по x
%
T=1;      % Заданный интервал времени
tau=h/3;  % Шаг по времени. 
          % Максимальное по модулю собственное значение матрицы А равно 3.
          % Поэтому:
          % Для устойчивости 'KIR' --- tau<=h/3  (нестрогий знак!)
          % Для 'Bim' --- tau<2*h/3 (строгий знак!)
%
method_type='Bim'; % Выбор метода. Допустимые значения:
                   % 1) 'KIR' --- Курант-Изаксон-Рис
                   % 2) 'Bim' --- Бим-Уорминг
%
u_init=[x+x.^3;    % Начальное условие
         x.^3;
        x.^2+1]; 
u_res=0*u_init;    % Итоговые значения
%
lambda=diag([0,2,3]); % Собственные значения матрицы А
U=[1 -1 -1;
   1 -2  1;
   1 -5  3];          % Базис из левых собственных векторов матрицы A
invU=[-1/4  2  -3/4;
      -1/2  1  -1/2;
      -3/4  1  -1/4]; % Обратная матрица к U
%
sigma=lambda*tau/h;
I=eye(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(method_type,'KIR'))
    u_work=U*u_init;
    for timesteps=1:T/tau
        for j=2:length(x)
            u_res(:,1)=u_work(:,1);
            u_res(:,j)=sigma*u_work(:,j-1)+(I-sigma)*u_work(:,j);
        end
        u_work=u_res;
    end
    u_res=invU*u_res;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(method_type,'Bim'))
    u_work=U*u_init;
    for timesteps=1:T/tau
        for j=3:length(x)
            u_res(:,1)=u_work(:,1);
            u_res(:,2)=sigma*u_work(:,1)+(I-sigma)*u_work(:,2);
            u_res(:,j)=(sigma-I)*sigma*u_work(:,j-2)/2-...
                             (sigma-2*I)*sigma*u_work(:,j-1)+...
                                  (I-(3*I-sigma)*sigma/2)*u_work(:,j);
        end
        u_work=u_res;
    end
    u_res=invU*u_res;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Отрисовка результатов:
%
zoom_x=[0,1];
%
figure(1)
hold on
plot(x,u_res(1,:),'b-');
plot(x,u_init(1,:),'r-');
xlabel('$x$','Interpreter','latex');
title('$u_1(x)$','Interpreter','latex');
xlim(zoom_x)
%
figure(2)
hold on
plot(x,u_res(2,:),'b-');
plot(x,u_init(2,:),'r-');
xlabel('$x$','Interpreter','latex');
title('$u_2(x)$','Interpreter','latex');
xlim(zoom_x)
%
figure(3)
hold on
plot(x,u_res(3,:),'b-');
plot(x,u_init(3,:),'r-');
xlabel('$x$','Interpreter','latex');
title('$u_3(x)$','Interpreter','latex');
xlim(zoom_x)