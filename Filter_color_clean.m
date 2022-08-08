%% KALMAN FILTER WITH COLOR V2

%% Build the single suspension
dt=0.001;% Timestep [s]
T=1000;% Total simulated time [s]
t=0:SYS.dt:T; % Time vector
N=numel(t); % number of iterations of the physical process

run singlesus_SS.m %Create the model

%Kalman_mode="Sequential"; % Sequential for sequential Kalman, anything else for the default mode
%Kalman_mode="Default"; % Sequential for sequential Kalman, anything else for the default mode
Kalman_mode="Steady-state"; %Steady-state Kalman filtering;


%% Get the inputs for the system

%Suspoint acceleration
as=0.000001*(.3*sin((sin(0.02*t+1).*1.*t)./((0.5*t-5).^2+1))+.8*cos(((sin(0.4*t+1)+1).*t./((t-200).^2+1))).*t./(t.^2+1)-5*exp(1./(((sin((sin(0.04*t+1)+1)*0.01.*t)).^2+3)))+3);
as=as-mean(as);
%M1 forces 
f1=0.001*(1*exp(-sin(.2*pi*2*t-300)).*exp(-0.1*sqrt((t-500).^2+0.1))-1*exp(-sin(2*pi*2*t-300)).*exp(-0.005*sqrt((t-15000).^2+0.1))); 

% Total input

u=[as; %Accelerations of the suspoint
   f1]; %Forces on M1

%Input as seen from the kalman state estimator, we don't know the HAM
%acceleration, so we assume it is zero plus some large error
uk=[zeros(1,N);
    f1];
   
%Process noise for the real system
w=transpose(mvnrnd(zeros(size(SYS.Q,1),1),SYS.Q,N));

%Measurement noise for the real system
v=transpose(mvnrnd(zeros(size(SYS.R,1),1),SYS.R,N));

%% PROPAGATE THE SYSTEM
%Init for all the arrays
K=zeros(SYS.size,2,N);
P=zeros(SYS.size,SYS.size,N);
X_pred=zeros(SYS.size,N);
X_minus=zeros(SYS.size,N);
X=zeros(SYS.size,N);
y=zeros(2,N);

% set the initial values
X(:,1)=10^-9*rand(SYS.size,1); %Some random initialization
P(:,:,1)=10^-15*eye(SYS.size);
y(:,1)=SYS.H*X(:,1)+v(:,1);

%Steady-state Kalman filtering
X_pred=zeros(SYS.size,N);

%Propagate the system forward in time
for i=1:(N-1)
    
    X(:,i+1)=SYS.F*X(:,i)+SYS.G*u(:,i)+w(:,i);
    y(:,i+1)=SYS.H*X(:,i+1)+v(:,i+1);
    
    %Propagate the predictions forward in time
    X_minus(:,i+1)=SYS.F*X_pred(:,i)+SYS.G*uk(:,i);
    P_minus=SYS.F*P(:,:,i)*transpose(SYS.F)+SYS.Qk; %Use the Kalman filter version of the process noise    
    
    switch(Kalman_mode)
        case("Sequential")
            if(i==1)
                disp("Sequential Kalman Filtering, interpret the Kalman gain as successive filters\n NOTE: assumes R is diagonal")
            end
%%%         % Sequential Kalman Filter
            X_pred(:,i+1)=X_minus(:,i+1);
            P(:,:,i+1)=P_minus;
            for j=1:size(y,1) %Iterate over every measurement
                Sigma=SYS.H(j,:)*P(:,:,i+1)*transpose(SYS.H(j,:))+SYS.R(j,j); % Sum of variances for weighted avg, should be a scalar
                if(Sigma<=0)
                    disp('Sigma<0, The variance is negative, aborting code')
                    return
                end
                K_seq=P(:,:,i+1)*transpose(SYS.H(j,:))/Sigma;

                %Perform measurement update
                X_pred(:,i+1)=X_pred(:,i+1)+K_seq*(y(j,i+1)-SYS.H(j,:)*X_pred(:,i+1));
                P(:,:,i+1)=(eye(SYS.size)-K_seq*SYS.H(j,:))*P(:,:,i+1)*transpose((eye(SYS.size)-K_seq*SYS.H(j,:)))+K_seq*SYS.R(j,j)*transpose(K_seq);
                K(:,j,i+1)=K_seq;
            end   
        case("Steady-state")
            if(i==1)
             disp("Only compute Steady-state filtering") 
            end
%%%         Steady-state Kalman filtering    
            X_pred(:,i+1)=(eye(SYS.size)-SYS.Kss*SYS.H)*(SYS.F*X_pred(:,i)+SYS.G*uk(:,i))...
                        +SYS.Kss*y(:,i+1);
            P(:,:,i+1)=(eye(SYS.size)-SYS.Kss*SYS.H)*P_minus*transpose((eye(SYS.size)-SYS.Kss*SYS.H))...
                        +SYS.Kss*SYS.R*transpose(SYS.Kss);
        otherwise
%%%         %Find the Kalman filter all at once
            Sigma=SYS.H*P_minus*transpose(SYS.H)+SYS.R; 
            K(:,:,i+1)=transpose(linsolve(Sigma,SYS.H*P_minus));

            %Measurement update
            X_pred(:,i+1)=X_minus(:,i+1)+K(:,:,i+1)*(y(:,i+1)-SYS.H*X_minus(:,i+1));
            P(:,:,i+1)=(eye(SYS.size)-K(:,:,i+1)*SYS.H)*P_minus*transpose((eye(SYS.size)-K(:,:,i+1)*SYS.H))...
                        +K(:,:,i+1)*SYS.R*transpose(K(:,:,i+1));
    end

end

