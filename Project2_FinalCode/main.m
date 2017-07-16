
clear all;
close all;
set(0,'defaulttextinterpreter','latex')

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


%%System parameters

A=[0 1 0;0 -0.87 43.22; 0 0.99 -1.34];
B=[0 0; -17.25 -1.58;-0.17 -0.25];
n=size(A,1);
m=size(B,2);
R1 = 0.3*eye(n);
R2= 0.01*eye(m);
Rq2=0.1*eye(n);
Rq3 = 0.05*eye(n);
Rq3(1,1)=0;
sigma = 0.5;
N=500;  % Number of trajectories for learning 

%% Learning the 'A', matrix


for itr=1:N
mdlA = sde(@DriftRateA, @DiffusionRate, 'StartState',[1 ; 0.5; -1]);    % Simulation model for dX = F(t,X)dt + G(t,X)dW

[X , time , Z] = simByEuler(mdlA, 100, 'DeltaTime', 0.01); % Simulation Solver

T=size(time,1);
%u=[sin(time),cos(time)];;
for i=1:T-1
     y(i,:)=X(i+1,:)-X(i,:);
 end
 
 for i=1:T-1
      phiA(i,:)= [X(i,:)]*(time(i+1)-time(i));
 end 
 
 for j=1:n
     W(:,:,j)=zeros(T-1,T-1);
     for i=1:T-1
         if X(i,j)~= 0
             W(i,i,j)=1/X(i,j)^2;
         end
     end
 end
 
 for j=1:n
     thetaA(j,:,itr)=(inv(phiA'*W(:,:,j)*phiA)*phiA'*W(:,:,j)*y(:,j))';
 end
end
Ahat = mean(thetaA,3);   %Estimation of A


%% Learning the 'B' matrix



for itr=1:N
mdlB = sde(@DriftRateB, @DiffusionRate, 'StartState',[1 ; 0.5; -1]);    % Simulation model for dX = F(t,X)dt + G(t,X)dW

[X , time , Z] = simByEuler(mdlB, 100, 'DeltaTime', 0.01); % Simulation Solver

 T=size(time,1);
 u=[sin(time),cos(time)];
for i=1:T-1
     yB(i,:)=X(i+1,:)-X(i,:)- X(i,:)*Ahat'*(time(i+1)-time(i));
 end
 
 for i=1:T-1
      phiB(i,:)= [u(i,:)]*(time(i+1)-time(i));
 end 
 
 
  for j=1:n
     W(:,:,j)=zeros(T-1,T-1);
     for i=1:T-1
         if X(i,j)~= 0
             W(i,i,j)=1/X(i,j)^2;
         end
     end
 end
 
 for j=1:n
     thetaB(j,:,itr)=(inv(phiB'*W(:,:,j)*phiB)*phiB'*W(:,:,j)*yB(:,j))';
 end
end
Bhat = mean(thetaB,3); %Estimation of B


%%Learning the diffusion coefficient 'sigma'


for itr=1:N
mdlSg = sde(@DriftRateSg, @DiffusionRate, 'StartState',[1 ; 0.5; -1]);    % Simulation model for dX = F(t,X)dt + G(t,X)dW

[X , time , Z] = simByEuler(mdlSg, 100, 'DeltaTime', 0.01); % Simulation Solver

 T=size(time,1);
 u=[sin(time),cos(time)];
 
 for i=1:T-1
     ySg(i,:)=X(i+1,:)-X(i,:);
 end
 for i=1:T-1
      xHat(i,:)= [X(i,:), u(i,:)]*(time(i+1)-time(i));
 end 
  
  for j=1:n
      tmp=0;
      for i=1:T-1
          tmp=tmp+(ySg(i,j)-[Ahat(j,:),Bhat(j,:)]*xHat(i,:)')^2/((time(i+1)-time(i))*X(i,j)^2);
      end
      sg(j)=tmp;
  end
 sgHat(itr) = sqrt(mean(sg)/(T-1));
end
sigmaHat=mean(sgHat); %Estimate of sigma


%Utilize the learned model dynamics for the Feedback control 

Bbar = Bhat*sqrt(inv(R2));
S=Bbar*Bbar';
Abar= Ahat+(1/2)*sigmaHat^2*eye(n);
Co=ctrb(Abar,Bbar);
unco=length(Abar)-rank(Co);
Cbar=sqrt(R1);
Ob = obsv(Abar,Cbar);
unob = length(Abar)-rank(Ob);
if (unco ~= 0) || (unob ~= 0)
   return;
end
[P,L,G] = care(Abar,Bbar,R1);

A2=Ahat+(3/2)*sigmaHat^2*eye(n)-S*P;
M2 = lyap(A2,Rq2);

A3=Ahat+(5/2)*sigmaHat^2*eye(n)-S*P;
M3 = lyap(A3,Rq3);

mdl = sde(@DriftRate, @DiffusionRate, 'StartState',[1 ; -0.5; -1]);    % Simulation model for dX = F(t,X)dt + G(t,X)dW
[X , time , Z] = simByEuler(mdl, 15000, 'DeltaTime', 0.001); % Simulation Solver


 for i=1:size(time,1)
 ctrl(i,:)= -inv(R2)*Bhat'*(P+(X(i,:)*M2*X(i,:)')*M2+ (X(i,:)*M3*X(i,:)')^2*M3)*X(i,:)';    
 end



%% Plot results

   figure(1);
   subplot(2,1,1)
   hold on
   plot(time', X);  
   title('State v. Time','fontsize',16); 
    xlabel('$t$, sec','fontsize',12)
    ylabel('$x(t)$','fontsize',12)
   legend('x_1(t)','x_2(t)','x_3(t)')
   hold off;
    
   subplot(2,1,2);hold on
   plot(time', ctrl);  
   title('Control v. Time','fontsize',16); 
    xlabel('$t$, sec','fontsize',12)
    ylabel('$u(t)$','fontsize',12)
    legend('\phi_1(x)','\phi_2(x)')
