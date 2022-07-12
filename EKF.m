function [mus,Sigs,Knorms_OSEM,Knorms_GS13] = EKF(U,Y,mu0,Sig0,f,g,A,C,Q,R) % just input system? Note must pass these as function handles
    n = numel(mu0);
    [p,T] = size(Y); % NOTE** Y should have shape T, and first entry goes unused
    
    mus = zeros(n,T);
    Sigs = zeros(n,n,T);
    mus_pred = zeros(n,T);
    Sigs_pred = zeros(n,n,T);
    Knorms_OSEM = zeros(1,T);
    Knorms_GS13 = zeros(1,T);
    
    mu = mu0;
    Sig = Sig0;
    
    for tt=1:(T-1) % Note different indexing
        % Store current belief
        mus(:,tt) = mu;
        Sigs(:,:,tt) = Sig;
        
        % Compute current dynamics Jacobian
        A_t = A(mu,U(:,tt));
        
        % Predict step
        mu_bar = f(mu,U(:,tt));
        Sig_bar = A_t*Sig*A_t' + Q;
        
        % Store predict beliefs
        mus_pred(:,tt) = mu_bar;
        Sigs_pred(:,:,tt) = Sig_bar;
        
        % Compute current measurement Jacobian
        C_t = C(mu_bar,U(:,tt));
        
        % Update step
        innov = Y(:,tt+1) - g(mu_bar,U(:,tt)); % Ensure angular outputs are wrapped
        
        K_t = (Sig_bar*C_t'*inv(C_t*Sig_bar*C_t' + R));
        Knorms_OSEM(tt+1) = norm(K_t(:,1:2),2);
        Knorms_GS13(tt+1) = norm(K_t(:,3),2); % Will generalize to higher dim GS13 outputs
        mu = mu_bar + K_t*innov;
        Sig = (eye(n) - K_t*C_t)*Sig_bar; 
    end
    
    % Store final belief and return
    mus(:,tt) = mu;
    Sigs(:,:,tt) = Sig;
    mus_pred(:,tt) = mu_bar; % These are indexed improperly
    Sigs_pred(:,:,tt) = Sig_bar;
    
    
    
end
        