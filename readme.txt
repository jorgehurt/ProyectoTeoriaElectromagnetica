%Universidad del Valle de Guatemala
%Jorge Hurtado, Fernando Sandoval, Mario Soto
%Proyecto de Simulacion
%Soluciones de Ecuacion de LaPlace
%Teoria Electromagnetica 1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% LaPlace en Coordenadas Rectangulares %%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problema 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIMER PROBLEMA INCISO A
1;
function [x,y,z,q,t,phi,efield,p] = poten(a,b,n)
  xp = -a:0.1:a;
  yp = -b:0.1:b;
  [x, y] = meshgrid(xp, yp);
  k = 0;
  Eo = 8.84.*10.^-12;
  z = x;
  for i = 1:n
    h = i*pi/b;
    fun = @(y)(2./(b.*sinh(h.*a))).*atan(y./a).*sin(h.*y);
    Fn = (integral(fun,0,b));
    Pot= Fn.* sinh(h.*x).*sin(h.*y);
    k += Pot;
  endfor
  phi = k;
  efield = gradient(-phi);
  [q,t] = gradient(-phi);
  p = Eo.*divergence(q,t);
endfunction 
a=3;
b=3;
Eo = 8.84.*10.^-12;
% potenciales
figure('Name','Grafica de Potencial de LaPlace en Rectangulares')
title('Grafica de Potencial de LaPlace en Rectangulares');
[x,y,z,q,t,phi,efield,p] = poten(a,b,2);
surf(x,y,phi);
colorbar
xlabel('x'), ylabel('y'), zlabel('Potencial')
grid on

% campos electricos 
figure('Name','Grafica de Campo Electrico de LaPlace en Rectangulares')
title('Grafica de Campo Electrico de LaPlace en Rectangulares');
quiver(q,t);
%quiver3(x,y,z,q,t,efield);
colorbar
xlabel('x'), ylabel('y'), zlabel('Campo Electrico')
grid on
% densidad de carga superficial
figure('Name','Densidad de Carga')
title('Densidad de Carga');
plot(y(1:61),(gradient(phi))(1:61))
xlim([-1 3])
xlabel('Y')
ylabel('Densidad de Carga')
grid on



%% PRIMER PROBLEMA INCISO B
1;
function [x,y,z,q,t,phi,efield,p] = poten(a,b,n)
  xp = -a:0.1:a;
  yp = -b:0.1:b;
  [x, y] = meshgrid(xp, yp);
  k = 0;
  Eo = 8.84.*10.^-12;
  e = 2.71828;
  z = x;
  for i = 1:n
    h = i*pi/b;
    j = i*pi*((a/b) + 1);
    d = i*pi;
    g = i*pi*(a/b);
    fun = @(y)((2.*y.^3)+5).*sin(h.*y);
    Fn = (integral(fun,0,b));
    An = ((e.^g+(-e.^-d))./(b.*sinh(j.*a))).*Fn;
    Bn = ((e.^d+(-e.^-g))./(b.*sinh(j.*a))).*Fn;
    Pot= ((An.*e.^((h.*x)) + Bn.*e.^((-h.*x)./b)).*sin(h.*y));
    k += Pot;
  endfor
  phi = k;
  efield = gradient(-phi);
  [q,t] = gradient(-phi);
  p = Eo.*divergence(q,t);
endfunction 
a=1;
b=1;
figure('Name','Grafica de Potencial de LaPlace en Rectangulares')
title('Grafica de Potencial de LaPlace en Rectangulares');
[x,y,z,q,t,phi,efield,p] = poten(a,b,2);
surf(x,y,phi);
colorbar
xlabel('x'), ylabel('y'), zlabel('Potencial')
grid on

% campos electricos 
figure('Name','Grafica de Campo Electrico de LaPlace en Rectangulares')
title('Grafica de Campo Electrico de LaPlace en Rectangulares');
%quiver(q,t);
quiver3(x,y,z,q,t,efield);
colorbar
xlabel('x'), ylabel('y'), zlabel('Campo Electrico')
grid on
% densidad de carga superficial
figure('Name','Densidad de Carga')
title('Densidad de Carga');
plot(y(1:21),(gradient(phi))(1:21))
xlim([-1 1])
xlabel('Y')
ylabel('Densidad de Carga')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problema 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SEGUNDO PROBLEMA INCISO A
1;
function [x,y,z,q,t,phi,efield,p] = poten(a,b,n)
  xp = -a:0.1:a;
  yp = -b:0.1:b;
  [x, y] = meshgrid(xp, yp);
  Eo = 8.84.*10.^-12;
  z=x;
  k=0;
  for i = 1:n
    h = i*pi/b;
    d = i*pi;
    g = i*pi*(a/b);
    fun = @(y) (2./b.*exp(-g)).*(atan(y./a)).*sin((h.*y));
    F = (integral(fun,0,b));
    Pot= F.*exp(-h.*x).*sin(h.*y);
    k += Pot;
  endfor
  phi = k;
  efield = gradient(-phi);
  [q,t] = gradient(-phi);
  p = Eo.*divergence(q,t);
endfunction 
a=1;
b=1;
figure('Name','Grafica de Potencial de LaPlace en Rectangulares')
title('Grafica de Potencial de LaPlace en Rectangulares');
[x,y,z,q,t,phi,efield,p] = poten(a,b,1);
surf(x,y,phi);
colorbar
xlabel('x'), ylabel('y'), zlabel('Potencial')
grid on
% campos electricos 
figure('Name','Grafica de Campo Electrico de LaPlace en Rectangulares')
title('Grafica de Campo Electrico de LaPlace en Rectangulares');
%quiver(q,t);
quiver3(x,y,z,q,t,efield);
colorbar
xlabel('x'), ylabel('y'), zlabel('Campo Electrico')
grid on
% densidad de carga superficial
figure('Name','Densidad de Carga')
title('Densidad de Carga');
plot(y(1:21),(gradient(phi))(1:21))
xlim([-1 1])
xlabel('Y')
ylabel('Densidad de Carga')
grid on

%% SEGUNDA PROBLEMA, INCISO B
1;
function [x,y,z,q,t,phi,efield,p] = poten(a,b,n)
  xp = -a:0.1:a;
  yp = -b:0.1:b;
  [x, y] = meshgrid(xp, yp);
  k = 0;
  Eo = 8.84.*10.^-12;
  e = 2.71828;
  z = x;
  for i = 1:n
    h = i*pi/b;
    j = i*pi.*((a/b) + 1);
    d = i*pi;
    g = i*pi*(a/b);
    fun = @(y)((2.*e.^(-g))./b).*((2.*y.^3)+5).*sin(h.*y);
    Fn = (integral(fun,0,b));
    Pot= Fn.* e.^(-h.*x).*sin(h.*y);
    k += Pot;
  endfor
  phi = k;
  efield = gradient(-phi);
  [q,t] = gradient(-phi);
  p = Eo.*divergence(q,t);
endfunction 
a=1;
b=1;
figure('Name','Grafica de Potencial de LaPlace en Rectangulares')
title('Grafica de Potencial de LaPlace en Rectangulares');
[x,y,z,q,t,phi,efield,p] = poten(a,b,20);
surf(x,y,phi);
colorbar
xlabel('x'), ylabel('y'), zlabel('Potencial')
grid on

% campos electricos 
figure('Name','Grafica de Campo Electrico de LaPlace en Rectangulares')
title('Grafica de Campo Electrico de LaPlace en Rectangulares');
%quiver(q,t);
quiver3(x,y,z,q,t,efield);
colorbar
xlabel('x'), ylabel('y'), zlabel('Campo Electrico')
grid on

% densidad de carga superficial
figure('Name','Densidad de Carga')
title('Densidad de Carga');
plot(y(1:21),(gradient(phi))(1:21))
xlim([-1 1])
xlabel('Y')
ylabel('Densidad de Carga')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% LaPlace en Coordenadas Esfericas %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problema 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grafica de Potencial de LaPlace en Coordenadas Esfericas
%Declaración de variables
R = 5;
e0=1;
 
%Creacion de Grafica
[s,theta] = meshgrid(0:0.1:10,0:pi./40:2.*pi);
[s,phie] = meshgrid(0:0.1:10,0:pi./80:pi);
PotIn =((3./(10.*e0)).*s.*cos(theta))+((-12./(35.*e0.*R.*R)).*((1./2).*(5.*(cos(theta).^3)-3.*cos(theta))));
PotOut = (((3.*R.*R.*R)./(10.*e0.*s.*s)).*cos(theta)+((-12.*(R^7))./(35.*e0.*R.*R.*(s.^4))).*((1./2).*(5.*(cos(theta).^3)-3.*cos(theta))));
PotGraf=PotIn.*(s<R)+PotOut.*(s>=R);
C = s.*theta;
%%[x,y] = pol2cart(theta,s);
[x,y,z]= sph2cart (phie,theta, s);
E=-gradient(PotGraf);
%Ploteo de Grafica en esfericas
figure('Name','Graficas de Potencial Esfericas de LaPlace en Esfericas Problema 1')
title('Grafica de Potencial Esfericas de LaPlace en Esfericas');
tab1 = ('Potencial' );
ax= axes('title','Potencial');
%%Grafica 1
subplot(2,2,1)
surf(-s,theta,PotGraf);
colorbar
xlabel('s'), ylabel('theta'), zlabel('Potencial')
grid on
%%Grafica 2
subplot(2,2,3)
contour(s,theta,PotGraf);
colorbar
xlabel('s'), ylabel('theta Potencial CN'), zlabel('Diferencial de Potencial')
grid on
%%Grafica 3
subplot(1,2,2)
[ds,da]=gradient(-PotGraf);
contour(s,theta, E);
hold on
colorbar
quiver3(s,theta,phie,ds,da,E);
xlabel('s'), ylabel('theta Campo Electrico'), zlabel('phie')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Problema 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grafica de Potencial de LaPlace en Coordenadas Esfericas
%Declaración de variables
R = 1;
e0=1;
P0=-5;
a=3;
b=6;
k=2; 
%Creacion de Grafica
[s,theta] = meshgrid(0:0.1:10,0:pi./40:2.*pi);
[s,phie] = meshgrid(0:0.1:10,0:pi./80:pi);

Pb=((-b.*k)./(e0)).*(1-(a./s))+((2.*k)./(15.*b.*e0)).*((s.*s)-((a.^5)./(s.^3))).*((3.*(cos(theta).^2))-(1))+((P0.*a)./(s));
Pc=((-b.*k.*(b-a))./(e0.*s))+((2.*k.*((b.^5)-(a.^5)))./(15.*b.*e0)).*(3.*(cos(theta).^2)-1)+((P0.*a)./(s));
Pa=P0;
PotGraf=Pa.*(s<a)+Pb.*(s<=b && s>=a)+Pc.*(s>b);


C = s.*theta;
%%[x,y] = pol2cart(theta,s);
[x,y,z]= sph2cart (phie,theta, s);
E=-gradient(PotGraf);
%Ploteo de Grafica en esfericas
figure('Name','Graficas de Potencial Esfericas de LaPlace en Esfericas Problema 1')
title('Grafica de Potencial Esfericas de LaPlace en Esfericas');
tab1 = ('Potencial' );
ax= axes('title','Potencial');
%%Grafica 1
subplot(2,2,1)
surf(s,theta,PotGraf);
colorbar
xlabel('s'), ylabel('theta'), zlabel('Potencial')
grid on
%%Grafica 2
subplot(2,2,3)
contour(s,theta,PotGraf);
colorbar
xlabel('s'), ylabel('theta Potencial CN'), zlabel('Diferencial de Potencial')
grid on
%%Grafica 3
subplot(1,2,2)
[ds,da]=gradient(-PotGraf);
contour(s,theta, E);
hold on
colorbar
quiver3(s,theta,phie,ds,da,E);
xlabel('s'), ylabel('theta Campo Electrico'), zlabel('phie')
hold off
%Grafica 4
%Densidad de carga superficial
%nuevo Plot de densidad
figure('Name','Densidad de Carga')
title('Densidad de Carga');
plot(y(1:101),(gradient(PotGraf))(1:101))
%xlim([-1 3])
xlabel('Y')
ylabel('Densidad de Carga')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% LaPlace en Coordenadas Cilindricas %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Grafica de Potencial de LaPlace en Coordenadas Cilindricas
%Declaración de variables
R = 2;
e0 = 1;
%Llamada a la grafica general
figure('Name','Grafica de Potencial de LaPlace en Cilindricas')
title('Grafica de Potencial de LaPlace en Cilindricas');
%Creacion de Grafica
[s,theta] = meshgrid(2:0.1:10,0:pi./50:2.*pi);
Z = (-(-e0).*s.*cos(theta).*(((R.*R)./(s.^2))-1));
C = s.*theta;
[x,y] = pol2cart(theta,s);
%Ploteo de Grafica
surf(x,y,Z,C);
colorbar
xlabel('x'), ylabel('y'), zlabel('potencial')
grid on