function f=lorenz_ext(~,X)
% Values of parameters
SIGMA = 10; R = 28; BETA = 8/3;
x=X(1); y=X(2); z=X(3);
Y= [X(4), X(7), X(10);
 X(5), X(8), X(11);
 X(6), X(9), X(12)];
f=zeros(9,1);
%Lorenz equation
f(1)=SIGMA*(y-x);
f(2)=-x*z+R*x-y;
f(3)=x*y-BETA*z;
%Linearized system
Jac=[-SIGMA, SIGMA, 0;
 R-z, -1, -x;
 y, x, -BETA];
%Variational equation
f(4:12)=Jac*Y;