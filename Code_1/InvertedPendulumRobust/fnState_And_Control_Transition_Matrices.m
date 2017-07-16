function [A,B] = fnState_And_Control_Transition_Matrices(x,u,du,dt)

global m;
global l;
global I;
global b;
global gr;


x1 = x(1,1);
x2 = x(2,1);

u1 = u(1,1);

A = zeros(2,2);


B = zeros(2,1);

A(1,1) = 0;
A(1,2)  = 1;
A(2,1)= -(gr/l)*cos(x1);
A(2,2)= -(b/I);


B(1,1)=0;
B(2,1)=1;



