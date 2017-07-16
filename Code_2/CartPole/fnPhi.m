function [xOut] = fnPhi(x,u,Horizon)
global mc;
global mp;
global l;
global gr;

for k = 1:(Horizon-1)
    
den = mc+mp*sin(x(3,k))^2;
xDDot = (1/den)*(u(k)+ mp*sin(x(3,k))*(l*x(4,k)^2 + gr*cos(x(3,k))));
      
xOut(k,1) = xDDot;
xOut(k,2) = xDDot*sin(x(3,k))^2 - gr*sin(x(3,k))*cos(x(3,k));
xOut(k,3) = -x(4,k)^2*sin(x(3,k));
end