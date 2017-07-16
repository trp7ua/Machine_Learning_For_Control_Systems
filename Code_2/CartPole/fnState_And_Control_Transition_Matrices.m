function [A,B] = fnState_And_Control_Transition_Matrices(x,u,du,dt)

global mcEst;
global mpEst;
global lEst;
global gr;

x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);
x4 = x(4,1);

den = mcEst+ (mpEst*sin(x3)^2);



A = zeros(4,4);


B = zeros(4,1);

A(1,2) = 1;
A(3,4)  = 1;


         
A(2,3) = (1/den^2)*((den*mpEst*(lEst*x4^2*cos(x3)+ gr*cos(2*x3)))-(u+mpEst*sin(x3)*(lEst*x4^2 + gr*cos(x3)))*(mpEst*sin(2*x3)));
A(2,4) = (1/den)*(2*mpEst*lEst*sin(x3)*x4);
A(4,3) = (1/(lEst*den)^2)*((lEst*den)*(u*sin(x3)-mpEst*lEst*x4^2*cos(2*x3)- (mcEst+mpEst)*gr*cos(x3)) + (u*cos(x3)+mpEst*lEst*x4^2*cos(x3)*sin(x3)+(mcEst+mpEst)*gr*sin(x3))*lEst*mpEst*sin(2*x3)); 
A(4,4) = -(1/den)*(mpEst*x4*sin(2*x3));
 
B(2,1) = 1/den;
B(4,1) = -(1/den)*(cos(x3)/lEst);





