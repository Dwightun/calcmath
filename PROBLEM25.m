function PROBLEM25
%
N=100;     % ����� ����� ����� �� x
h=1/N;    % ��� ����� �� x
x=linspace(0,1,N+1); % ���� ����� �� x
%
T=1;      % �������� �������� �������
tau=h/3;  % ��� �� �������. 
          % ������������ �� ������ ����������� �������� ������� � ����� 3.
          % �������:
          % ��� ������������ 'KIR' --- tau<=h/3  (��������� ����!)
          % ��� 'Bim' --- tau<2*h/3 (������� ����!)
%
method_type='Bim'; % ����� ������. ���������� ��������:
                   % 1) 'KIR' --- ������-�������-���
                   % 2) 'Bim' --- ���-�������
%
u_init=[x+x.^3;    % ��������� �������
         x.^3;
        x.^2+1]; 
u_res=0*u_init;    % �������� ��������
%
lambda=diag([0,2,3]); % ����������� �������� ������� �
U=[1 -1 -1;
   1 -2  1;
   1 -5  3];          % ����� �� ����� ����������� �������� ������� A
invU=[-1/4  2  -3/4;
      -1/2  1  -1/2;
      -3/4  1  -1/4]; % �������� ������� � U
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
% ��������� �����������:
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