function [x] = fnsimulate(xo,u_new,Horizon,sigma,dt)

global m;
global l;
global I;
global b;
global gr;

x = xo;
for k = 1:(Horizon-1)
          
      Fx(1,1) = x(2,k);
      Fx(2,1) = -(b/I)*x(2,k)-(gr/l)*sin(x(1,k));
      G_x(1,1) = 0 ;
      G_x(2,1) = 1;
x(:,k+1) = x(:,k) + Fx * dt + G_x * u_new(:,k) * dt + sqrt(sigma)*randn*dt;
end