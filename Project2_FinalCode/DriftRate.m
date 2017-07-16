function F = DriftRate(t,x)
global A;
global B;
global R1;
global R2;
global P;
global M2;
global M3;
global sigma;
global Ahat;
global Bhat;
global sigmaHat;

u=-inv(R2)*Bhat'*(P+(x'*M2*x)*M2+ (x'*M3*x)^2*M3)*x;
F=A*x+B*u;
