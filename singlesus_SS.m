%% Routine to make a state space for a single suspension for Kalman filtering
%
SYS.note={'State space matrices for the augmented model for Kalman filtering, all units are SI',...
    'Inputs are suspoint acceleration and top mass forces',...
    'process noise inputs are encoded in W, measurement noise inputs are incoded in V',...
    'The A,B,C,D,Qc,Rc,Wc,Vc matrices represent the continuous system.',...
    'The F,G,H,Q,R,W,V represent the discrete system, discretized by hand with a zero-order hold',...
    'The dynamical system is modeled up to suspoint velocity'};

%% Variable initializations
%%% Simulation variables %%%
f_lowpass=0.001; %Hz, all integrals are just low pass filters in disguise

%% Define the dynamical system

%%% Physical variables %%%
g=9.81; % gravitational acceleration [m/s^2]
L=5; %length of the pendlum stage [m]
m=1; %Mass of the pendulum [kg]
b=0.1; %Viscous damping coefficient [Ns/m]

if exist('var_as','var') == 0
    var_as=1; %Variance of the acceleration noise [m/s^2]^2
end

%      xs    x1    xs'    x1' 
DYN.A=[0,    0,    1,      0;   % xs'
       0,    0,    0,      1;   % x1'
       0,    0,    0,      0;   % xs''
      g/L, -g/L,   b/m,   -b/m]; % x1''

%      as   f1
DYN.B=[0,    0;   % xs'
       0,    0;   % x1'
       1,    0;   % xs''
       0,   1/m]; % xs''
   
%       xs   x1    xs'   x1'    
DYN.C=[ 0,    0,    1,    0; % y1 (GS13 velocity measurement)
        -1,   1,    0,    0];% y2 (OSEM velocity measurement)

    
%           ws
DYN.Wc=[    0; % xs'
            0;   % x1'
            0;   % xs''
            0]; % x1''
        
%
DYN.Vc= [];

%Continuous time covariances
DYN.Qc=DYN.Wc*transpose(DYN.Wc);
DYN.Rc=DYN.Vc*transpose(DYN.Vc);

% get the size of the state space
DYN.size=size(DYN.A,1);

DYN=Discretize(DYN,dt,f_lowpass);

%%% The Kalman filter model assumes high covariance for the acceleration of
%%% the table
DYN.Qk=DYN.Q;
DYN.Qk(3,3)=var_as*dt;

DYN.Qck=DYN.Qc;
DYN.Qck(3,3)=var_as;
          
%% Create the noise states
%The noise states are handled by subroutines that can be tested and
%debugged independently

GS13=define_GS13_states(dt,f_lowpass);
OSEM=define_OSEM_states(dt,f_lowpass);

%% Append the states to create a bigger state space model

% System ins and outs
SYS.A=blkdiag(DYN.A,GS13.A,OSEM.A);
SYS.B=[ DYN.B;
        GS13.B;
        OSEM.B];
SYS.C=[DYN.C,GS13.C,OSEM.C];

% Noise inputs
SYS.Wc=blkdiag(DYN.Wc,GS13.Wc,OSEM.Wc);
SYS.Vc=[DYN.Vc,GS13.Vc,OSEM.Vc];

%Covariances for process and measurement noises
SYS.Qc=SYS.Wc*transpose(SYS.Wc);
SYS.Rc=SYS.Vc*transpose(SYS.Vc);

%Discretize the system
SYS.size=size(SYS.A,1);
SYS=Discretize(SYS,dt,f_lowpass);

%Add the covariance for the 'infinite' process noise used for the kalman
%filter
SYS.Qk=SYS.Q;
SYS.Qk(1:DYN.size,1:DYN.size)=DYN.Qk;

%Add the covariance for the 'infinite' process noise used for the kalman
%filter
SYS.Qck=SYS.Qc;
SYS.Qck(1:DYN.size,1:DYN.size)=DYN.Qck;

% Pass other important information about the model
SYS.dt=dt; 
SYS.f_lowpass=f_lowpass;

%Define state-spaces for both the Kalman filter and the input measurements
[SYS.Pc,SYS.Kc_t]=icare(SYS.A,transpose(SYS.C),SYS.Qck,SYS.Rc); %Get the steady-state Kalman gain

%% Define the GS13 variables
function GS13=define_GS13_states(dt,f_lowpass)
%%% FLAGS FOR DEBUGGING %%%
debug_flag=0;
    % Definitions for pink noise
     f1=1; f2=f1/4; f3=f2/4; f4=f3/4; f5=f4/4;
     g1=1; g2=2*g1; g3=2*g2; g4=2*g3; g5=2*g4;
     gtot=(9.5*10^-10)/sqrt(2);

     f_pink=[f1;f2;f3;f4;f5];
     g_pink=gtot*[g1;g2;g3;g4;g5];
     n_pink=numel(f_pink);
     
    % Definitions for white (measurement) noise
     g_white=2.2*(10^-11)/sqrt(2);
     
    %              nu_pink        nu_pink_int     nu_pink_int_int
     GS13.A=[-2*pi*diag(f_pink), zeros(n_pink,1), zeros(n_pink,1);  % nu_pink'
               ones(1,n_pink)  ,      0         ,        0       ;  % nu_pink_int'
               zeros(1,n_pink) ,      1         ,        0       ]; % nu_pink_int_int'

     %              as              f1    
     GS13.B=[zeros(n_pink,1), zeros(n_pink,1);  % nu_pink'
                    0       ,        0       ;  % nu_pink_int'
                    0       ,        0       ]; % nu_pink_int_int' 

     %              nu_pink        nu_pink_int     nu_pink_int_int
     GS13.C=[zeros(1,n_pink),           0       ,          1      ; % y1 (GS13 velocity measurement)
             zeros(1,n_pink),           0       ,          0      ];% y2 (OSEM velocity measurement)

     %             wpink
     GS13.Wc=[2*pi*f_pink.*g_pink; % nu_pink'
                    0           ; % nu_pink_int'
                    0           ];% nu_pink_int_int'

     %          nu_white GS13
     GS13.Vc=[     g_white        ; % y1 (GS13 velocity measurement)
                     0           ];% y2 (OSEM velocity measurement)
                 
    %Continuous time covariances
    GS13.Qc=GS13.Wc*transpose(GS13.Wc);
    GS13.Rc=GS13.Vc*transpose(GS13.Vc);

    %%% Define some Quality-of-life outputs
     GS13.freq=logspace(-2,2,1000); 
     GS13.intended_noise = SEI_sensor_noise('GS13meas',GS13.freq); 
     GS13.size=size(GS13.A,1);
    
    %%% Discretize the equations %%%
    GS13=Discretize(GS13,dt,f_lowpass);
    
    %%% DEBUGGING AREA %%%
    if(debug_flag)
        n_debug=0;
        for i=1:n_pink
            n_debug=n_debug+g_pink(i)./((1i)*freq./f_pink(i)+1);
        end
        resp_debug=sqrt((abs(n_debug)./(2*pi*freq).^2).^2+GS13.Vc(1).^2);
        
        %simulate the noise
        T=1000;
        t=0:dt:T;
        N=numel(t);
        X=zeros(GS13.size,N);
        y=zeros(2,N);
        w=transpose(mvnrnd(zeros(size(GS13.Q,1),1),GS13.Q,N));
        v=transpose(mvnrnd(zeros(size(GS13.R,1),1),GS13.R,N));
        disp(size(w))
        for i=1:(N-1)
            X(:,i+1)=GS13.F*X(:,i)+w(:,i);
            y(:,i+1)=GS13.H*X(:,i+1)+v(:,i+1);
        end
        % Calculate the ASD
        [asd_noise,freq2]=asd2(y(1,:), dt,  5, 3, @hann, 'symmetric');
        
        %%% PLOT %%%
        figure,test=loglog(freq2,asd_noise,freq,GS13.intended_noise.*(2*pi*freq),...
                            freq,sqrt(2)*resp_debug); %Don't forget sqrt(2) for the ASD
        grid on; title("SS debugging for GS13 Noise"); ylabel("Amplitude [m/s\cdot\surd{Hz}]");
        xlabel("Frequency [Hz]"); set(test(2:3),"linewidth",3); set(gca,"fontsize",14);
        xlim([10^-2,10^2])
        legend("Simulated Noise", "GS13 Noise", "Fitted dynamical model");
        figure,plot(t,y(1,:))
    end
end
 %% Define the OSEM variables
function OSEM=define_OSEM_states(dt,f_lowpass)
%%% FLAGS FOR DEBUGGING %%%
debug_flag=0;
    % Definitions for pink noise
     f1=12.5; f2=f1/4; f3=f2/4;
     g1=2/2.1; g2=2*g1; g3=2*g2;
     gtot=(1.7*10^-11)/sqrt(2);

     f_pink=[f1;f2;f3];
     g_pink=gtot*[g1;g2;g3];
     n_pink=numel(f_pink);
    
    % definitions for brown noise
     g_brown=4.62*10^-10/sqrt(2);
    
    % definitions for white (measurement) noise
     g_white=2.82*10^-11/sqrt(2);
    %              nu_pink        nu_brown     
     OSEM.A=[-2*pi*diag(f_pink), zeros(n_pink,1);  % nu_pink'
               zeros(1,n_pink),        0        ]; % nu_brown'

     %              as              f1    
     OSEM.B=[zeros(n_pink,1), zeros(n_pink,1);  % nu_pink'
                    0       ,        0       ];  % nu_brown'

     %            nu_pink           nu_brown     
     OSEM.C=[zeros(1,n_pink),           0       ; % y1 (GS13 velocity measurement)
             ones(1,n_pink),            1       ];% y2 (OSEM velocity measurement)

     %             wpink
     OSEM.Wc=[2*pi*f_pink.*g_pink, zeros(n_pink,1); % nu_pink'
                    0           ,    g_brown     ]; % nu_pink_int'


     %          nu_white GS13
     OSEM.Vc=[        0           ; % y1 (GS13 velocity measurement)
                  g_white        ];% y2 (OSEM velocity measurement)

    %Continuous time covariances
    OSEM.Qc=OSEM.Wc*transpose(OSEM.Wc);
    OSEM.Rc=OSEM.Vc*transpose(OSEM.Vc);
    
    %%% Define some Quality-of-life outputs
     OSEM.freq=logspace(-2,2,1000);
     OSEM.intended_noise = HAM_SUS_OSEM_noise('L',OSEM.freq); 
     OSEM.size=size(OSEM.A,1); 
    
    %%% Discretize the equations %%%
    OSEM=Discretize(OSEM,dt,f_lowpass);
    
    %%% DEBUGGING AREA %%%
    if(debug_flag)
        n_debug=0;
        for i=1:n_pink
            n_debug=n_debug+g_pink(i)./((1i)*OSEM.freq./f_pink(i)+1);
        end
        resp_debug=sqrt(abs(n_debug).^2+OSEM.Wc(end,end)^2./(2*pi*OSEM.freq).^2+OSEM.Vc(end,end).^2); 
        
        %%% simulate the noise %%%
        T=1000;
        t=0:dt:T;
        N=numel(t);
        X=zeros(OSEM.size,N);
        y=zeros(2,N);
        w=transpose(mvnrnd(zeros(size(OSEM.Q,1),1),OSEM.Q,N));
        v=transpose(mvnrnd(zeros(size(OSEM.R,1),1),OSEM.R,N));
        for i=1:(N-1)
            X(:,i+1)=OSEM.F*X(:,i)+w(:,i);
            y(:,i+1)=OSEM.H*X(:,i+1)+v(:,i+1);
        end
        % Calculate the ASD
        [asd_noise,freq2]=asd2(y(2,:), dt,  5, 3, @hann, 'symmetric'); %The OSEM noise affects y2
        
        %%% PLOT %%%
        figure,test=loglog(freq2,asd_noise,OSEM.freq,OSEM.intended_noise,...
                                 OSEM.freq,sqrt(2)*resp_debug); %Don't forget sqrt(2) for the ASD
        grid on; title("SS debugging for OSEM Noise"); ylabel("Amplitude [m/\surd{Hz}]");
        xlabel("Frequency [Hz]"); set(test(2:3),"linewidth",3); set(gca,"fontsize",14);
        xlim([10^-2,10^2])
        legend("Simulated Noise", "OSEM Noise", "Fitted dynamical model"); 
        figure,plot(t,y(2,:))
    end
end
%% Routine to discretize equations 
function Disc=Discretize(Cont,dt,f_lowpass)
% Converts the continuous-time equations into discrete time, with an
% 'exact' approximation. 
%Differential equations get integrated, 
% white noise scales as 1/sqrt(dt), which affects process and measurement noises
% differently (one gets integrated and the other does not)

% Copy the array
Disc=Cont;

% All integrals are lowpass filters in disguise
Disc.A=Disc.A-2*pi*f_lowpass*eye(Cont.size); 

% Discretize the dynamics
Disc.F=expm(Disc.A*dt);
Disc.H=Cont.C;
Mat=Disc.A^0*dt;
for i=1:100
    Mat=Mat+(-Disc.A)^i*dt^(i+1)/factorial(i+1);
end
Disc.G=Disc.F*(Mat)*Cont.B;

% Discretize the noises
Disc.W=Cont.Wc*sqrt(dt); % The process noise gets integrated
Disc.V=Cont.Vc/sqrt(dt); % The measurement noise doesn't get integrated

%Discrete time covariances
Disc.Q=Disc.Qc*dt;
Disc.R=Disc.Rc/dt;
end