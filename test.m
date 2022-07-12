close all
clear all
format compact
clc
%% Set constants
% Add random seed*
As = .1;
fss = linspace(.2,10,50);
%% Plot sensor noise ASDs
% figure()
% n_OSEM = HAM_SUS_OSEM_noise('L', fss);
% n_GS13 = SEI_sensor_noise('GS13meas', fss);
% plot(fss,n_OSEM)
% hold on
% plot(fss,n_GS13)
% xlabel('Frequency')
% ylabel('Spectral Noise Density (m/rtHz)')
% legend('OSEM','GS13')
% title('Comparison of Sensor Noise Spectra')
% saveas(gcf,'sens_noise.png')
% hold off

%% Execute ideal, frequency-informed filter
% Assumes only suspension point motion
% Assumes suspension point motion is sinusoidal
% Executes ideal extended kalman filter (EKF) once per every input frequency
% Collects euclidean norm of errors and plot

for i=1:numel(fss)
    fs = fss(i);
    % Simulate motion
    sus = doubsus(As,fs); % Implement brownian process noise
    [eqRed_1,eqRed_2] = sus.EOMs;
    [sols,xs] = sus.solveEOMs(eqRed_1,eqRed_2); % Defaults to zero ICs
    t = sols.x;
    thetas_thetaDots = sols.y;
    % sus.plotsols(sols,xs)
    % model.plotsols(sols,xs,true,true)
    
    Xc = sus.theta2disp(thetas_thetaDots,xs); % Convert theta-state and suspension point movement to displacement-state
    % Xc = [xs, x1, x2]
    [ncont,T] = size(Xc);
    U = zeros(1,T);

    % Take measurements
    % fs = sus.fs;
    L_1 = sus.L_1;
    L_2 = sus.L_2;
    GS13 = sens('GS13meas');
    OSEM = sens('OSEM');
    Y_GS13 = GS13.measureAtFrequency(fs,Xc,U);
    Y_OSEM = OSEM.measureAtFrequency(fs,Xc,U);

    % FILTER for one specific frequency
    % X_k+1 = A_k*X_k + B_k*U_k + W_k
    % Y_k = C_k*X_k + V_k
    % Cov(W_k) = Q
    % Cov(V_k) = R
    
    ndisc = 2*ncont;
    Y_combo = [Y_OSEM; Y_GS13];

    f_sus = @(Xt,Ut) sus.flin(Xt,Ut,fs);
    A_sus = @(Xt,Ut) sus.Ajaclin(Xt,Ut,fs);
    Q_sus = eye(ndisc);

    g_OSEM = @(Xt,Ut) OSEM.g(Xt,Ut); % Assumes xs = Ut
    C_OSEM = @(Xt,Ut) OSEM.Cjaclin(Xt,Ut);
%     R_OSEM = OSEM.informedNoise(fs,2);
    g_GS13 = @(Xt,Ut) GS13.g(Xt,Ut);
    C_GS13 = @(Xt,Ut) GS13.Cjaclin(Xt,Ut);
%     R_GS13 = GS13.informedNoise(fs,1);

    g_combo = @(Xt,Ut) [g_OSEM(Xt,Ut);g_GS13(Xt,Ut)];
    C_combo = @(Xt,Ut) [C_OSEM(Xt,Ut);C_GS13(Xt,Ut)];
%     R_combo = [R_OSEM, zeros(2,1);
%          zeros(1,2), R_GS13];
    R_combo = r*eye(3);

    mu0 = zeros(ndisc,1); % Add some disturbance**
    Sig0 = eye(ndisc); % State space dimension
    [mus,Sigs,Knorms_OSEM,Knorms_GS13] = EKF(U,Y_combo,mu0,Sig0,f_sus,g_combo,A_sus,C_combo,Q_sus,R_combo);

    % Analyze results
    x0 = Xc(1,:);
    x1 = Xc(2,:);
    x2 = Xc(3,:);

    mu_x0 = mus(1,:);
    mu_x1 = mus(2,:);
    mu_x2 = mus(3,:);

    x0_meas = Y_GS13;
    x1rel_meas = Y_OSEM(1,:);
    x2rel_meas = Y_OSEM(2,:);
    
    % Compute norm of residual errors
    % Compare results from simply using OSEMs to using filtered OSEM+GS13
    x1_OSEM_resid(i) = norm(x1rel_meas-x1,2);
    x2_OSEM_resid(i) = norm(x2rel_meas-x2,2);
    x1_EKF_resid(i) = norm(mu_x1-x1,2);
    x2_EKF_resid(i) = norm(mu_x2-x2,2);
    
    

%     figure()
%     plot(t,x1rel_meas)
%     hold on
%     plot(t,x2rel_meas)
%     plot(t,Y_GS13)
%     legend('OSEM_1','OSEM_2','GS13')
%     xlabel('Time')
%     ylabel('Y_t')
%     title('Measurements')
%     hold off
% 
%     figure()
%     plot(t,Knorms_OSEM)
%     hold on
%     plot(t,Knorms_GS13)
%     legend('EKF OSEM weighting','EKF GS13 weighting')
%     xlabel('Time')
%     ylabel('||K_t||_2')
%     hold off
% 
%     figure()
%     plot(t, x1)
%     hold on
%     plot(t, x2)
%     plot(t,x1rel_meas)
%     plot(t,x2rel_meas)
%     plot(t,mu_x1)
%     plot(t,mu_x2)
%     legend('x_1 truth','x_2 truth','x_1 OSEM estimate','x_2 OSEM estimate','x_1 EKF estimate','x_2 EKF estimate')
%     xlabel('Time')
%     ylabel('State')
%     hold off
end

% Plot euclidean norm of errors
figure()
plot(fss,x1_OSEM_resid)
hold on
plot(fss,x2_OSEM_resid)
plot(fss,x1_EKF_resid)
plot(fss,x2_EKF_resid)
legend('x1 OSEM residual','x2 OSEM residual','x1 EKF residual','x2 EKF residual')
xlabel('Frequency')
ylabel('Residual (m)')
title('Residual error of estimates at different frequencies, uniform noise baseline')
saveas(gcf,'ResidErr_same.png')
hold off
