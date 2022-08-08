%% KALMAN FILTER WITH COLOR V2

%% Build the single suspension

dt=1/4096;%timestep [s]
T=2000;% Total simulated time [s]
t=0:dt:(T-dt); % Time vector
N=numel(t); % number of iterations of the physical process

run singlesus_SS.m %Create the actual model
%% Judiciously create the state-space models for running a simulation

%Real system, all states are outputs. The noises are appended (but they
%                     inputs , noises
Contss.real=ss(SYS.A,[ SYS.B , SYS.Wc ],eye(SYS.size),[]); %We need to add the noises as inputs to the noise states.

%Construct the state estimator
%                                                  inputs    measurements
Contss.kalman=ss(SYS.A-transpose(SYS.Kc_t)*SYS.C,[ SYS.B ,transpose(SYS.Kc_t)],eye(SYS.size),[]); %We add the measurements as inputs to the state estimator

%Construct the equivalent discrete-time object for filtering
Discss.real=c2d(ss(SYS.A,SYS.B,SYS.C,[]),dt);
[P,K]=idare(transpose(Discss.real.A),transpose(SYS.H),SYS.Qk,SYS.R);
Discss.kalman=ss(Discss.real.A*(eye(SYS.size)-transpose(K)*SYS.H),...
    [Discss.real.B, Discss.real.A*transpose(K)],...
    eye(SYS.size),[],dt);
%% Define the inputs for the state-space system
%%% REAL INPUTS %%%
%Suspoint acceleration
[~, ~, ~, as] = GenerateTimeSeries(@HAM_asd_generator, dt, T);
%M1 forces 
f1=0.0001*(1*exp(-sin(.2*pi*2*t-300)).*exp(-0.1*sqrt((t-500).^2+0.1))-1*exp(-sin(2*pi*2*t-300)).*exp(-0.005*sqrt((t-15000).^2+0.1))); 
% Total input

u=[transpose(as); %Accelerations of the suspoint
       f1]; %Forces on M1

%Input as seen from the kalman state estimator, we don't know the HAM
%acceleration, so we assume it is zero
uk=[zeros(1,N);
    f1];
   
%%% NOISE INPUTS %%%
% (White) noise scales as 1/sqrt(dt)


%Process noise, input for the real system. Using Wc implies we use unit
%variance noise
w=transpose(mvnrnd(zeros(size(SYS.Wc,2),1),eye(size(SYS.Wc,2)),N))/sqrt(dt);

%Measurement noise for the real system
v=transpose(mvnrnd(zeros(size(SYS.Rc,2),1),SYS.Rc,N))/sqrt(dt);
%% Simulate the system
X0=10^-9*rand(SYS.size,1); %Some random initialization

X=transpose(lsim(Contss.real,[u;w],t,X0));     % System simulation
y=SYS.C*X+v;                  % Simulate the GS13 and OSEM measurements
X_pred=transpose(lsim(Contss.kalman,[uk;y],t));% Kalman filter the measurements
X_pred_disc=transpose(lsim(Discss.kalman,[uk;y],t)); %Discrete Kalman filter
%% CALCULATE ASD FOR RECORD KEEPING

%%% Frequency domain %%%
begin=floor(N/3); %only get the tail of the behavior

% SUSPOINT POSITION
[ASD.Xs,freq2]=asd2(X(1,begin:end), dt,  9, 3, @hann, 'symmetric'); %true motion
[ASD.Xs_pred,~]=asd2(X_pred(1,begin:end), dt,  9, 3, @hann, 'symmetric');%filter prediction
[ASD.Xs_res,~]=asd2(X(1,begin:end)-X_pred(1,begin:end), dt,  9, 3, @hann, 'symmetric'); %residual
[ASD.Xs_disc_res,~]=asd2(X(1,begin:end)-X_pred_disc(1,begin:end), dt,  9, 3, @hann, 'symmetric'); %discrete residual

% SUSPOINT VELOCITY
[ASD.Vs,freq2]=asd2(X(3,begin:end), dt,  9, 3, @hann, 'symmetric'); %true motion
[ASD.Vs_pred,~]=asd2(X_pred(3,begin:end), dt,  9, 3, @hann, 'symmetric');%filter prediction
[ASD.Vs_res,~]=asd2(X(3,begin:end)-X_pred(3,begin:end), dt,  9, 3, @hann, 'symmetric'); %residual
[ASD.Vs_disc_res,~]=asd2(X(3,begin:end)-X_pred_disc(3,begin:end), dt,  9, 3, @hann, 'symmetric'); %discrete residual

% X1
[ASD.X1,freq2]=asd2(X(2,begin:end), dt,  9, 3, @hann, 'symmetric'); %true motion
[ASD.X1_pred,~]=asd2(X_pred(2,begin:end), dt,  9, 3, @hann, 'symmetric');%filter prediction
[ASD.X1_res,~]=asd2(X(2,begin:end)-X_pred(2,begin:end), dt,  9, 3, @hann, 'symmetric'); %residual
[ASD.X1_OSEM,~]=asd2(y(2,begin:end), dt,  9, 3, @hann, 'symmetric'); %relative motion
[ASD.X1_OSEM_res,~]=asd2(y(2,begin:end)-X(2,begin:end), dt,  9, 3, @hann, 'symmetric'); %relative motion_residual
[ASD.X1_disc_res,~]=asd2(X(2,begin:end)-X_pred_disc(2,begin:end), dt,  9, 3, @hann, 'symmetric'); %residual

ASD.freq2=freq2;

%Response to suspoint motion.
%NOTE: this is done by hand, needs to be automated.
Xs_to_X1=squeeze(freqresp(Contss.real(2,1),2*pi*GS13.freq)).*transpose(-(2*pi*GS13.freq).^2);


%% PLOT RESULTS

% SUSPOINT POSITION
figure, comparison=loglog(freq2,ASD.Xs,freq2,ASD.Xs_res,freq2,ASD.Xs_disc_res,...
                    GS13.freq,GS13.intended_noise,'--',OSEM.freq,OSEM.intended_noise,'--');
                
grid on; xlim([10^-2,10^2]); ylim([10^-15,10^-3])
set(gca,'fontsize',15)
legend('Real Xs motion','Estimator residual','Estimator residual (discrete)','GS13 Noise','OSEM Noise')
xlabel('Frequency [Hz]')
ylabel('Amplitude [m/\surd{Hz}]')
title('Suspoint position comparison')

% SUSPOINT VELOCITY
figure, comparison=loglog(freq2,ASD.Vs,freq2,ASD.Vs_res,freq2,ASD.Vs_disc_res,...
                    GS13.freq,GS13.intended_noise.*(2*pi.*GS13.freq),'--');               
grid on; xlim([10^-2,10^2]); ylim([10^-15,10^-3])
set(gca,'fontsize',15)
legend('Real suspoint velocity','Estimator residual','Estimator residual (discrete)','GS13 Velocity Noise')
xlabel('Frequency [Hz]')
ylabel('Amplitude [m/(s\cdot\surd{Hz})]')
title('Suspoint velocity comparison')

% X1 POSITION
figure, comparison=loglog(freq2,ASD.X1,freq2,ASD.X1_res,freq2,ASD.X1_disc_res,...
                    freq2,ASD.X1_OSEM_res,GS13.freq,GS13.intended_noise,OSEM.freq,OSEM.intended_noise,GS13.freq,...
                    transpose(abs(Xs_to_X1)).*GS13.intended_noise,'--');
grid on; xlim([10^-2,10^2])
set(gca,'fontsize',15)
legend('Real X1 motion', 'Estimator residual','Estimator residual (discrete)','OSEM_residual','GS13 Noise','OSEM Noise','GS13 feedforward noise')
xlabel('Frequency [Hz]')
ylabel('Amplitude [m/\surd{Hz}]')
title('M1 position comparison')

%%% TIME SERIES %%%

% % X1 motion
 figure, mass_displacement=plot(t/60,X_pred(2,:)-X(2,:),t/60,X_pred_disc(2,:)-X(2,:));
 grid on, xlim([0,T/60]);
 set(gca,'fontsize',15)
 legend('Kalman prediction residual','Kalman prediction (discrete) residual')
 xlabel('Time [min]')
 ylabel('Amplitude [m]')

% Suspoint Velocity
figure, suspoint_displacement=plot(t/60,X_pred(3,:)-X(3,:),t/60,X_pred_disc(3,:)-X(3,:),t/60,y(1,:)-X(3,:),'--');
grid on, xlim([0,T/60]);
set(gca,'fontsize',15)
legend('Kalman prediction','Kalman prediction (discrete)','GS13 measurement')
xlabel('Time [min]')
ylabel('Amplitude [m]')
