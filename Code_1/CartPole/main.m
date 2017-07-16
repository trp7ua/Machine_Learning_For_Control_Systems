%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Iterative Optimal Control for Cart Pole Dynamics  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Machine Learning  for Control Systems    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE8803 Fall 2015                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Tanmay Rajpurohit              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
set(0,'defaulttextinterpreter','latex')

global mc;
global mp;
global l;
global gr;

% Cart Pole Model Parameter
% masses in Kgr
 mc = 1;
 mp = 0.01;

 %Gravitation Acceleration
 gr = 9.81;
 
% length parameters in meters
 l = 0.25;
 
% Horizon 
Horizon = 300; % 1.5sec
% Number of Iterations
num_iter = 100

% Discretization
dt = 0.01;

% Weight in Final State:
Q_f = zeros(4,4);
Q_f(2,2) = 5;
Q_f(3,3) = 5;
Q_f(4,4) = 5;

% Weight in the Control:
R = 0.1 ;

% Initial Configuration:
xo = zeros(4,1);

% Initial Control:
u_k = zeros(1,Horizon-1);
du_k = zeros(1,Horizon-1);


% Initial trajectory:
x_traj = zeros(4,Horizon);
 

% Target: 
p_target(1,1) = 2;
p_target(2,1) = 0;
p_target(3,1) = pi;
p_target(4,1) = 0;


% Learning Rate:
gamma = 0.08
 
 
for k = 1:num_iter

%------------------------------------------------> Linearization of the dynamics
%------------------------------------------------> Quadratic Approximations of the cost function 
for  j = 1:(Horizon-1)
    
     [l0,l_x,l_xx,l_u,l_uu,l_ux] = fnCost(x_traj(:,j), u_k(:,j), j,R,dt);
      q0(j) = dt * l0;
      q_k(:,j) = dt * l_x;
      Q_k(:,:,j) = dt * l_xx;
      r_k(:,j) = dt * l_u;
      R_k(:,:,j) = dt * l_uu;
      P_k(:,:,j) = dt * l_ux; 
    
    [dfx,dfu] = fnState_And_Control_Transition_Matrices(x_traj(:,j),u_k(:,j),du_k(:,j),dt);
   
    A(:,:,j) = eye(4,4) + dfx * dt;
    B(:,:,j) = dfu * dt;  
end

%------------------------------------------------> Find the controls
Vxx(:,:,Horizon)= Q_f;
Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target); 
V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target); 


%------------------------------------------------> Backpropagation of the Value Function
for j = (Horizon-1):-1:1
     
   H = R_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * B(:,:,j);
   G = P_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j);  
   g = r_k(:,j) +  B(:,:,j)' * Vx(:,j+1);
   
   
   inv_H = inv(H);
   L_k(:,:,j)= - inv_H * G;
   l_k (:,j) = - inv_H *g;  
   

   Vxx(:,:,j) = Q_k(:,:,j)+ A(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j) + L_k(:,:,j)' * H * L_k(:,:,j) + L_k(:,:,j)' * G + G' * L_k(:,:,j);
   Vx(:,j)= q_k(:,j) +  A(:,:,j)' *  Vx(:,j+1) + L_k(:,:,j)' * g + G' * l_k(:,j) + L_k(:,:,j)'*H * l_k(:,j);
   V(:,j) = q0(j) + V(j+1)   +   0.5 *  l_k (:,j)' * H * l_k (:,j) + l_k (:,j)' * g;
end 


%----------------------------------------------> Find the controls
dx = zeros(4,1);
for i=1:(Horizon-1)    
   du = l_k(:,i) + L_k(:,:,i) * dx;
   dx = A(:,:,i) * dx + B(:,:,i) * du;  
   u_new(:,i) = u_k(:,i) + gamma * du;
end

u_k = u_new;


%---------------------------------------------> Simulation of the Nonlinear System
[x_traj] = fnsimulate(xo,u_new,Horizon,dt);
[Cost(:,k)] =  fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
x1(k,:) = x_traj(1,:);
 

fprintf('Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
 
 
end



   time(1)=0;
   for i= 2:Horizon
    time(i) =time(i-1) + dt;  
   end

      

   figure(1);
   subplot(3,2,1)
   hold on
   plot(time,x_traj(1,:),'linewidth',4);  
   plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
   title('$x$','fontsize',20); 
    xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
   
    subplot(3,2,2);hold on;
   plot(time,x_traj(2,:),'linewidth',4); 
   plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
   title('$\dot{x}$','fontsize',20);
    xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
   
   

    
    subplot(3,2,3);hold on
   plot(time,x_traj(3,:),'linewidth',4); 
   plot(time,p_target(3,1)*ones(1,Horizon),'red','linewidth',4)
   title('$\theta$','fontsize',20)
     xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
   subplot(3,2,4);hold on
   plot(time,x_traj(4,:),'linewidth',4); 
   plot(time,p_target(4,1)*ones(1,Horizon),'red','linewidth',4)
   title('$\dot{\theta}$','fontsize',20)
     xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
    subplot(3,2,5);hold on
   plot(Cost,'linewidth',2); 
    xlabel('Iterations','fontsize',20)
   title('Cost','fontsize',20);
   %save('DDP_Data');