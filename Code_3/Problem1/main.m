%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  %%                                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Machine Learning  for Control Systems    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE8803 Fall 2015                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Tanmay Rajpurohit              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
set(0,'defaulttextinterpreter','latex')

%------Parameters of the model-----

A= 1.5;
B=0.1;
sigmaEp= 0.01;
Q=10;
R=1;
x0 = 2;     %--Intial State
N=251;           %--- Time horizon
theta= -1;       %-- Intial value of parameter
trajN=250;       %--- Number of trajectories for computation of expectation
gamma=0.01;     %-- Learning rate


itr=0;          %Iteration counter 

gradient=100;   %-- Some intial value to perform first iteration
while (abs(gradient) >  1 )  %-- convergence criterion
itr=itr+1;

% ---- Computation of trajectories
for j=1:trajN
    x(1,j)=x0;
    u(1,j)=theta*x0;
    for k=1:N-1
        x(k+1,j)=A*x(k,j)+B*u(k,j)+sigmaEp*normrnd(0,1);
        u(k+1,j)=theta*x(k+1,j);
    end
end

%---- Computation of Reward and storing as function of iteration
for j=1:trajN
    reward(j)=x(:,j)'*x(:,j)*Q+u(:,j)'*R*u(:,j);
end
J(itr)=mean(reward);

% ---- Computation of trajectories after perturbation in theta
delTheta=0.01;
tmpTheta = theta+delTheta;
for j=1:trajN
    x(1,j)=x0;
    u(1,j)=tmpTheta*x0;
    for k=1:N-1
        x(k+1,j)=A*x(k,j)+B*u(k,j)+sigmaEp*normrnd(0,1);
        u(k+1,j)=tmpTheta*x(k+1,j);
    end
end

 %---- Computation of Reward after perturbation
for j=1:trajN
    reward(j)=x(:,j)'*x(:,j)*Q+u(:,j)'*R*u(:,j);
end
delJ = J(itr)-mean(reward);

%-- Computing gradient
gradient=delJ/delTheta;

%---Updating the policy parameters
if abs(gradient) > 10000
    theta = theta+ sign(gradient);
elseif abs(gradient) > 1000
    theta = theta+ 0.5*sign(gradient);
elseif abs(gradient) > 100
    theta = theta+0.1*sign(gradient);
else
    theta=theta+gamma*gradient;
end

%-- Storing the value of parameters as function of iteration
KFD(itr)=theta;
end

%--MATLAB global optimal LQR gains
[K,S,E]=dlqr(A,B,Q,R);

%%--- Ploting the Results
figure(1)
hold on;
plot([1:itr],KFD);
plot(-K*ones(1,itr),'red','linewidth',4)
title('Plot of $K_{FD}$ versus Iterations','fontsize',16)
xlabel('Iteration','fontsize',12)
legend('K_{FD}','K_{LQR}')
hold off;

figure(2)
plot(J);
title('Reward versus Iterations','fontsize',16)
xlabel('Iteration','fontsize',12)

