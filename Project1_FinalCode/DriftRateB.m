function F = DriftRate(t,x)
global A;
global B;

u=[sin(t);cos(t)];
F=A*x+B*u;
