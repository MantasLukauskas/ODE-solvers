function f=rab_fab_ext(~,X)
% Values of parameters
gama = 0.1;
alfa = 0.2715;
SIGMA = 10; R = 28; BETA = 8/3;
x=X(1); y=X(2); z=X(3);
Y= [X(4), X(7), X(10);
 X(5), X(8), X(11);
 X(6), X(9), X(12)];
f=zeros(9,1);
%Lorenz equation
f(1)= y*(z-1+x*x) + gama*x;
f(2)= x*(3*z+1-x*x) + gama*y;
f(3)= -2*z*(alfa + x*y);
%Linearized system
Jac=[2*x*y+gama, z-1+x*x, y;
 3*z+1-3*x*x, gama, 3*x;
 -2*z*y, -2*z*x, -2*alfa];
%Variational equation
f(4:12)=Jac*Y;