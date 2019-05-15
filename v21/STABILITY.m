function STABILITY
%
% Выбор метода:
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
x=linspace(-12,1,101);
y=linspace(-7,7,101);
R=zeros(length(x),length(y));
[X,Y]=ndgrid(x,y);
gamma=1+sqrt(2)/2;
for i=1:length(x)
   for j=1:length(y)
          z=x(i)+1i*y(j);
          if(strcmp(method_type,'MidPoint'))
             R(i,j)=abs((1+z/2)/(1-z/2));
          end
          if(strcmp(method_type,'HH'))
             R(i,j)=abs((1+z/2+z^2/12)/(1-z/2+z^2/12));
          end
          if(strcmp(method_type,'KB'))
             R(i,j)=abs((1+z/2+z^2/10+z^3/120)/(1-z/2+z^2/10-z^3/120));
          end
          if(strcmp(method_type,'SDIRC'))
             R(i,j)=abs((1+(1-2*gamma)*z+0.5*(1-4*gamma+2*gamma^2)*z^2)/((1-gamma*z)^2));
          end
   end
end
figure(1)
hold on
xlabel('$\mathrm{Real}(z)$','Interpreter','latex')
ylabel('$\mathrm{Imag}(z)$','Interpreter','latex')
title('$|R(z)|$','Interpreter','latex')
contourf(X,Y,R,linspace(0,1,11),'ShowText','on')