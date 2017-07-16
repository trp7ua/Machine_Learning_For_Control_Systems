function [x] = fnsimulateEst(xo,u_new,Horizon,dt)

global mcEst;
global mpEst;
global lEst;
global gr;


x = xo;

for k = 1:(Horizon-1)
    

den = mcEst+mpEst*sin(x(3,k))^2;

      Fx(1,1) = x(2,k);
      Fx(2,1) = (1/den)*(mpEst*sin(x(3,k))*(lEst*x(4,k)^2 + gr*cos(x(3,k))));
      Fx(3,1) = x(4,k);
      Fx(4,1) = -(1/(lEst*den))*(mpEst*lEst*x(4,k)^2*cos(x(3,k))*sin(x(3,k)) + (mcEst+mpEst)*gr*sin(x(3,k)));
   
         
      G_x(1,1) = 0;
      G_x(2,1) = 1/den;
      G_x(3,1) = 0;
      G_x(4,1) = - cos(x(3,k))/(lEst*den);
      
x(:,k+1) = x(:,k) + Fx * dt + G_x * u_new(:,k) * dt;
end