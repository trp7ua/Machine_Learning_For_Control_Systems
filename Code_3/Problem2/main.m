%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
sigmaEp= 0.001;
sigma=0.001;
Q=0.1;
R=0.01;
x0 = 2;     %--Intial State
N=251;           %--- Time horizon
theta= -5;       %-- Intial value of parameter
trajN=500;       %--- Number of trajectories for computation of expectation
%gamma=0.01;     %-- Learning rate


itr=0;          %Iteration counter 

gradient=100;   %-- Some intial value to perform first iteration


while (abs(gradient) >  0.025 )  %-- convergence criterion
itr=itr+1;
    
% ---- Computation of trajectories
for j=1:trajN
    x(1,j)=x0;
    controlNoise(:,j)=normrnd(zeros(1,N),ones(1,N));
    u(1,j)=theta*x0+sigma*controlNoise(1,j);
    for k=1:N-1
        x(k+1,j)=A*x(k,j)+B*u(k,j)+sigmaEp*normrnd(0,1);
        u(k+1,j)=theta*x(k+1,j)+ sigma*controlNoise(k+1,j);
    end
end

%---- Computation of reward per trajectory and partial derivative of
%-----log-probability
for j=1:trajN
    reward(j)=x(:,j)'*Q*x(:,j)+u(:,j)'*R*u(:,j);
    gradProb(j)=(1/sigma)*(x(:,j)'*controlNoise(:,j));
end

%--- Computing Baseline and gradient
baseLine=mean(reward.*(gradProb.^2))/mean((gradProb.^2));
gradient = mean((reward - baseLine).*gradProb);


%---Updating the policy parameters
if abs(gradient) > 10
     theta = theta - sign(gradient);
elseif abs(gradient) > 1
    theta=theta - 0.5*sign(gradient);
else
    theta=theta - gradient;
end

%-- Storing the parameters and reward as function of iteration
KRPG(itr)=theta;
J(itr)=mean(reward);
end

%--MATLAB global optimal LQR gains
[K,S,E]=dlqr(A,B,Q,R);

%%--- Ploting the Results
figure(1)
hold on;
plot([1:itr],KRPG);
plot(-K*ones(1,itr),'red','linewidth',4)
title('Plot of $K_{RPG}$ versus Iterations','fontsize',16)
xlabel('Iteration','fontsize',12)
legend('K_{RPG}','K_{LQR}')
hold off;

figure(2)
plot(J);
title('Reward per Iterations','fontsize',16)
xlabel('Iteration','fontsize',12)

