%% Routine to make a quick Kalman filter

%% Variable initializations

%Physical variables
g=9.81; % gravitational acceleration [m/s^2]
L=5; %length of the pendlum stage [m]
m=1; %Mass of the pendulum [kg]
b=1; %Viscous damping coefficient [Ns/m]

%Simulation variables
dt=0.0001; % Timestep [s]
T=30;% Total simulated time [s]
t=0:dt:T;
N=numel(t);

% The state is defined as [xs ; x1; (x1)'] 
F=[ 1,        0,       0;
    0,        1,      dt;
    g*dt/L, -g*dt/L, 1-b*dt/m 
    ];
%The input is defined as [displacements in xs; forces in x1]
G=[ 1, 0;
    0,  0;
    0, dt/m
    ];
%The measurements are defined as [GS13; OSEM]
H=[  1 ,0, 0;
    -1, 1, 0
    ];
% Process noise covariance matrix, we assume that there is input noise in
% all hidden states. 
var_xs=10^10; %process covariance for the position x_s, assumed to be essentially infinite
var_x1=0;% process covariance for the position x_1, I don't know what to put in here
var_x1dot=0; %process covaricance for the velocity, this one is easily understood as random forces

Q=[var_xs,0,0;
    0,var_x1,0;
    0,0,var_x1dot
    ];
% Measurement noise covariance matrix
var_GS13=.0001; %Measurement covariance for the GS13 (assumed white here)
var_OSEM=.0002; %Measurement covariance for the OSEM (assumed white here)

R=[var_GS13, 0;
    0, var_OSEM
    ];
%% Creating random noise processes
nGS13=normrnd(0,sqrt(var_GS13),1,N);
nOSEM=normrnd(0,sqrt(var_OSEM),1,N);


%% Creating a fictitious input signal
xs_in=10*exp(-sqrt(t)-sin(5*t));
f1_in=1*exp(-sin(2*pi*2*t-3)).*exp(-sqrt((t-5).^2+0.1))-1*exp(-sin(2*pi*2*t-3)).*exp(-sqrt((t-15).^2+0.1)); 

xs_actual=0.1*exp(-cos(1.3*t-3)).*exp(-sqrt((t-5).^2+0.1));%Actual xs, that should be observed by the GS13

%Actual input notation
u_actual= zeros(2,N);
u_actual(1,1:(N-1))=xs_actual(2:N)-xs_actual(1:(N-1));
u_actual(2,:)=f1_in;

%Model input notation
u_model= zeros(2,N);
u_model(1,:)=0*xs_in;
u_model(2,:)=f1_in;

xGS13=xs_actual+nGS13-xs_actual(1); %Predefined measurement from the GS13

xOSEM(1)=0; % initialize the first OSEM measurement, 

%Measurement moving stuff
y=zeros(2,N);
y(1,:)=xGS13;
y(2,:)=xOSEM;

%% Initialize the Kalman filter predictions
K=zeros(3,2,N);
P=zeros(3,3,N);
P(:,:,1)=[10^10, 0, 0;
            0, 10^10,0;
            0, 0,10^10
            ];
X_pred=zeros(3,N);
X_pred2=zeros(3,N);
X_minus=zeros(3,N);
X_syst=zeros(3,N);
y_low=zeros(2,N);
tau=1/100; %in seconds
X_syst(:,1)=[0.001;0.1;0.2];
%% Iterating on the predictions
for i=1:(N-1)
    % Advance the actual system
    X_syst(:,i+1)=F*X_syst(:,i)+G*u_actual(:,i);
    % Measure the actual system (OSEM only)
    y(2,i+1)=H(2,:)*X_syst(:,i+1)+nOSEM(i+1);
    y(1,i+1)=H(1,:)*X_syst(:,i+1)+nGS13(i+1);
    % Advance the prediction
    X_minus(:,i+1)=F*X_pred(:,i)+G*u_model(:,i);
    % Find the new filter
    P_minus=F*P(:,:,i)*transpose(F)+Q;
    K(:,:,i+1)=P_minus*transpose(H)/(H*P_minus*transpose(H)+R);
    % Filter the prediction
    P(:,:,i+1)=(eye(3)-K(:,:,i+1)*H)*P_minus;
    %P(:,:,i+1)=(P(:,i+1)+transpose(P(:,i+1)))/2;
    X_pred(:,i+1)=(eye(3)-K(:,:,i+1)*H)*X_minus(:,i+1)+K(:,:,i+1)*y(:,i+1);
    X_pred2(:,i+1)=(eye(3)-K2*H)*(F*X_pred2(:,i)+G*u_model(:,i))+K2*y(:,i+1);
    y_low(:,i+1)=(1-dt/tau)*y_low(:,i)+dt/tau*y(:,i+1);
end
%% Plotting
K(:,:,i+1)
% Plot the predicted X1 motion
fin=floor(2*N/3);
figure, comparison=plot(t(1:fin),y(1,1:fin)+y(2,1:fin),t(1:fin),y_low(1,1:fin)+y_low(2,1:fin),t(1:fin),X_pred(2,1:fin),t(1:fin),X_pred2(2,1:fin),t(1:fin),X_syst(2,1:fin),'--');
set(comparison,'linewidth',2)
grid on
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('Unfiltered measurement', 'Single pole Lowpassed measurement (100 Hz)','Kalman Filtered Measurement','Kalman Filtered Measurement (steady state)','True Motion')
title('X1 motion')
set(gca,'fontsize',14)

figure, comparison3=plot(t(1:fin),y(1,1:fin),t(1:fin),y_low(1,1:fin),t(1:fin),X_pred(1,1:fin),t(1:fin),X_syst(1,1:fin));
set(comparison,'linewidth',2)
grid on
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('Unfiltered measurement', 'Single pole Lowpassed measurement (100 Hz)','Kalman Filtered Measurement','True Motion')
title('Xsuspoint motion')
set(gca,'fontsize',14)

figure, comparison2=plot(t(1:fin),y_low(1,1:fin)+y_low(2,1:fin)-X_syst(2,1:fin),t(1:fin),X_pred(2,1:fin)-X_syst(2,1:fin),t(1:fin),X_pred2(2,1:fin)-X_syst(2,1:fin));
set(comparison2,'linewidth',3)
grid on
xlabel('Time [s]')
ylabel('Amplitude [m]')
legend('Single pole Lowpassed measurement (100 Hz)','Kalman Filtered Measurement','Steady-State Kalman Filtered Measurement')
title('Residual comparison')
set(gca,'fontsize',14)