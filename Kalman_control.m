%% Control by using a state estimator
% This script is trying to use our dynamical system to test different
% control schemes with state estimators. It can be used to test the
% robustness of our Kalman control scheme.

%It works in conjunction with the Simulink diagram
%kalman_control_simulink.slx
%Ideally we populate all of the simulink diagram and use it to more
%efficiently implement the controls that we would like to test.
%% Define the dynamical system and its discretization parameters

dt=1/4096;%ADC sample time
T=2000;% Total simulated time [s]
t=0:dt:(T-dt); % Time vector
N=numel(t); % number of iterations of the physical process

%Create the dynamical model with [A,B,C,D] matrices
singlesus_SS

%% Populate the Simulink diagram

%% Dynamical system
dyn=ss(DYN.A+0*2*pi*eye(DYN.size)*f_lowpass,DYN.B,eye(DYN.size),[]); %Defined in continuous-space

%% GS13 noise model
%NOTE: In this case, different from the Kalman estimator version, we incorporate
%and feed-forward the white measurement noise to the output.

%                pink         white                              pink            white  
gs13=ss(GS13.A,[GS13.Wc,zeros(GS13.size,1)],GS13.C(1,:),[zeros(1,size(GS13.Wc,2)), 1  ]);

%% OSEM noise model
%NOTE: In this case, different from the Kalman estimator version, we incorporate
%and feed-forward the white measurement noise to the output

%              pink&brown      white                          pink&brown         white  
osem=ss(OSEM.A,[OSEM.Wc,zeros(OSEM.size,1)],OSEM.C(2,:),[zeros(1,size(OSEM.Wc,2)), 1  ]);

%% Measurement Model

%Matrix that combines the state space and the OSEM and GS13 noises
%       Output of Dyn   ,     add noise
C_meas=[    DYN.C       ,      eye(2)    ];

%% Estimator (pick one of the options below)
%ins: meas, drives
%outs:estim_meas, estim_states.

%%%KALMAN ESTIMATOR%%%
%Define the modeled system, which is seen by the Kalman estimator.
%NOTE Suspoint acceleration is not seen by the estimator
%                         inputs           noise
Modeled_ss=ss(SYS.A,[ SYS.B(:,2:end)   ,  SYS.Wck], SYS.C,[]);
kalm_disc=kalmd(Modeled_ss,eye(size(SYS.Wck,2)),SYS.Rc,dt);

%Grab only the extimated states
estimator=kalm_disc(size(SYS.C,1)+(1:size(DYN.A,1)),:);
%% Control related things:

%Signal Blending: We select signals from the ordered pool to lowpass and
%highpass
%              Sensors       Estimator
%            yGS13 yOSEM   xs  x1  xs' x1' 
% lowselector= [0,     1,    0,  0,  0,  0];%signal to feed lowpass
% highselector=[0,     0,    0,  1,  0,  0];%signal to feed highpass
% 
% 
% lp=ss(0); %Lowpass filter, defined directly
% hp=1-lp;%Highpass filter, this is supposed to be 1-lp
% 
% % controller
% % I am actually not very good at loop shaping, should I ask someone about
% % this?
% %zpk(0,2*pi*[1,1.2],-0)

% Use a linear-quadratic regulator (because it is easy to implement)

[~,K]=idare(DYN.F,DYN.B(:,2),1E5*[0,1,0,0]'*[0,1,0,0],1);
%              Sensors       Estimator
%            yGS13 yOSEM   xs  x1  xs' x1' 
lowselector= [0,     1,    0,  0,  0,  0];%signal to feed lowpass
%highselector=[0,     0,          K      ];%signal to feed highpass
highselector=[1,     0,    0,  0,  0,  0];%signal to feed highpass
lp=ss(0); %Lowpass filter, defined directly
hp=1-lp;%Highpass filter, this is supposed to be 1-lp

control=ss(-1);

% Create the State-Space (OPEN LOOP)
control_gain=0; %Open loop
SS.OL=linearize('kalman_control_simulink');

% Create the State-Space (CLOSED LOOP)
control_gain=1;
SS.CL=linearize('kalman_control_simulink');
%% SIMULATION

%INPUTS

%Suspoint acceleration
[~, ~, ~, as] = GenerateTimeSeries(@HAM_asd_generator, dt, T);
%M1 forces 
f1=0*10^-4*(1*exp(-sin(.2*pi*2*t-300)).*exp(-0.1*sqrt((t-500).^2+0.1))-1*exp(-sin(2*pi*2*t-300)).*exp(-0.005*sqrt((t-15000).^2+0.1))); 

%NOISE INPUTS
w=transpose(mvnrnd(zeros(size(SYS.Wc,2),1),eye(size(SYS.Wc,2)),N))/sqrt(dt);
v=transpose(mvnrnd(zeros(size(SYS.Rc,2),1),SYS.Rc,N))/sqrt(dt);

%TOTAL INPUTS
u=[as';f1;w(2,:);v(1,:);w(3:end,:);v(2,:)];

%SIMULATED OUTPUT
X.OL=lsim(SS.OL,u,t); %Open Loop
X.CL=lsim(SS.CL,u,t); %Closed Loop

[~,in,out]=In_Out_Parser(struct(),'kalman_control_simulink'); %Grab the channel names for easier access
%% GET ASD

%Open loop ASD
[ASD.OL.x1,freq2]=asd2(X.OL(:,out.m1.position), dt,  9, 3, @hann, 'symmetric');
[ASD.OL.v1,~]=asd2(X.OL(:,out.m1.velocity), dt,  9, 3, @hann, 'symmetric');

%Closed Loop ASD
[ASD.CL.x1,~]=asd2(X.CL(:,out.m1.position), dt,  9, 3, @hann, 'symmetric');
[ASD.CL.v1,~]=asd2(X.CL(:,out.m1.velocity), dt,  9, 3, @hann, 'symmetric');

%Estimator residual (Open loop)
[ASD.RES.x1,~]=asd2(X.OL(:,out.estim.m1.position)-X.OL(:,out.m1.position), dt,  9, 3, @hann, 'symmetric');
[ASD.RES.v1,~]=asd2(X.OL(:,out.estim.m1.velocity)-X.OL(:,out.m1.velocity), dt,  9, 3, @hann, 'symmetric');

%Input acceleration
[ASD.IN.as,~]=asd2(u(1,:), dt,  9, 3, @hann, 'symmetric');

%FAKE OSEM DAMPING
[ASD.RES.osem,~]=asd2(X.CL(:,out.OSEM.meas)-X.CL(:,out.m1.position), dt,  9, 3, @hann, 'symmetric');
%% PLOT
figure, figu=loglog(freq2,ASD.IN.as./(2*pi*freq2).^2,...
               freq2,ASD.OL.x1,freq2, ASD.RES.osem,'--',...
               freq2,ASD.CL.x1,freq2,ASD.RES.x1,'--',...
               OSEM.freq,OSEM.intended_noise,'--');
grid on; xlim([10^-2,10^2]); ylim([10^-15,10^-3])
set(gca,'fontsize',15)
set(figu,'linewidth',3)
legend('Suspension point displacement','Open-Loop X1 motion','OSEM Residual','Closed-Loop X1 motion','Estimator residual','OSEM Noise')
xlabel('Frequency [Hz]')
ylabel('Amplitude [m/\surd{Hz}]')
title('M1 position comparison')
%%