function STABILITY
%
% Выбор метода:
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
x=linspace(-12,1,101);
y=linspace(-7,7,101);
R=zeros(length(x),length(y));
[X,Y]=ndgrid(x,y);
gamma=1-sqrt(2)/2;
for i=1:length(x)
   for j=1:length(y)
          z=x(i)+1i*y(j);
          if(strcmp(method_type,'Euler'))
             R(i,j)=abs(1/(1-z));
          end
          if(strcmp(method_type,'Rado3'))
             R(i,j)=abs((1-z/3)/(1-2*z/3+z^2/6));
          end
          if(strcmp(method_type,'Rado5'))
             R(i,j)=abs((1+2*z/5+z^2/20)/(1-3*z/5+3*z^2/20-z^3/60));
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