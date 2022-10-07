%% CHANGELOG
    % REV0
    %% 
    % * Changed variable M_* to M_*
    % * Added x_0, M_0, u_*
    % * NOTE: x_0 hard-coded as a function of t
    %% 
    % REV1
    %% 
    % * Implemented optical control actuators dynamics
    %%
    % REV2
    %%
    % * Added plant model and connection to GS-13
    % REV3 (01JUN22)
    % * Converted to class ('doubsus.m')
    % * Pared down commentary
    % REV4 (TBD)
    % * Changed OSEM input model to be EM force from wall

classdef doubsus
    properties
        n = 2; % Dimension of base state** INCONSISTENT
        dt = 1/(64e3); % Sampling rate is ~64K *
        
        % Solver params
        X0 % Initial state
        trange
        
        
        % Seismic/suspension point input parameters
        As % amplitude of seismic disturbance
        fs % Freq of ""

        % Optical actuator input parameters (currently turned off)
        Ao1 = 0;
        fo1 = 0;
        Ao2 = 0;
        fo2 = 0;

        % Pendulum parameters
        g = 9.8;
        f0_1 % Natural frequency, see Clark 4.1
        f0_2 % Properties don't allow reference to other properties*
        L_1 % Properties don't allow reference to other properties* perhaps use a constructor to handle
        L_2
        mo_1 % mass of optic
        mt_1 % mass of remaining test mass
        mo_2 % mass of optic
        mt_2 % mass of remaining test mass
        M_1
        M_2
    end
    
    methods
        function obj = doubsus(As,fs) % Default constructor
            % Inputs are seismic/suspension point input parameters
            obj.As = As;
            obj.fs = fs;
%             obj.dt = .1; % Sampling time
            
            % Solver params
            obj.X0 = zeros(2*obj.n,1); % Initial state
            obj.trange = [0 50];


            % Optical actuator input parameters (currently turned off)
            obj.Ao1 = 0;
            obj.fo1 = 0;
            obj.Ao2 = 0;
            obj.fo2 = 0;

            % Pendulum parameters
            obj.g = 9.8;
            obj.f0_1 = .685; % Natural frequency, see Clark 4.1
            obj.f0_2 = obj.f0_1;
            obj.L_1 = obj.g/(obj.f0_1^2);
            obj.L_2 = obj.g/(obj.f0_2^2);
            obj.mo_1 = .5; % mass of optic
            obj.mt_1 = 1.5; % mass of remaining test mass
            obj.mo_2 = obj.mo_1;
            obj.mt_2 = obj.mt_1;
            obj.M_1 = obj.mo_1 + obj.mt_1;
            obj.M_2 = obj.mo_2 + obj.mt_2;
        end
        
        function r = test(obj)
            r=1;
        end
        
        function [eqRed_1,eqRed_2] = EOMs(obj)
            syms t theta_t_1(t) theta_t_2(t) %As fs Ao1 Ao2 fo1 fo2 L_1 L_2 mo_1 mo_2 mt_1 mt_2 g % x_0(t) --- not sure if t needs explicit definition theta_1(t) theta_2(t) u_1(t) u_2(t)
            
            As = obj.As;
            fs = obj.fs;
            Ao1 = obj.Ao1;
            fo1 = obj.fo1;
            Ao2 = obj.Ao2;
            fo2 = obj.fo2;
            g = obj.g;
%             f0_1 = obj.f0_1;
%             f0_2 = obj.f0_2;
            L_1 = obj.L_1;
            L_2 = obj.L_2;
            mo_1 = obj.mo_1;
            mt_1 = obj.mt_1;
            mo_2 = obj.mo_2;
            mt_2 = obj.mt_2;
            M_1 = obj.M_1;
            M_2 = obj.M_2;
            
            %% 
            % Define the displacements of the double pendulum in Cartesian coordinates.

            % Can implement seismic tilt disturbance here with an affine rotation*
            % x_0 = @(t) As*sin(2*pi*fs*t); % Fcn handles not diff()
            x_0 = As*sin(2*pi*fs*t);
            y_0 = 0;
            x_t_1 = x_0 + L_1*sin(theta_t_1); % The (t) after x_0 seems to not matter
            y_t_1 = y_0 + -L_1*cos(theta_t_1);
            x_t_2 = x_t_1 + L_2*sin(theta_t_2);
            y_t_2 = y_t_1 - L_2*cos(theta_t_2);

            % OSEM displacement/input
            % u_1 = @(t) Ao1*sin(2*pi*fo1*t);
            % u_2 = @(t) Ao2*sin(2*pi*fo2*t);
            u_1 = Ao1*(sin(2*pi*fo1*t)+1); % Offset s.t. only displacing mass in one direction
            u_2 = Ao2*(sin(2*pi*fo2*t)+1);
            %% 
            % Find the velocities by differentiating the displacements with respect to time 
            % using the <docid:symbolic_ug#btwol6i diff> function.

            vx_t_1 = diff(x_t_1);
            vy_t_1 = diff(y_t_1);
            vx_t_2 = diff(x_t_2);
            vy_t_2 = diff(y_t_2);
            %% 
            % Find the accelerations by differentiating the velocities with respect to time.

            ax_t_1 = diff(vx_t_1);
            ay_t_1 = diff(vy_t_1);
            ax_t_2 = diff(vx_t_2);
            ay_t_2 = diff(vy_t_2);
            %% Step 2: Define Equations of Motion
            % Define the equations of motion based on Newton's laws.
            % 
            % First, specify the tension of the first rod as $T_1$, and the tension of the 
            % second rod $T_2$.

            syms T_1 T_2
            %% 
            % Next, construct free-body diagrams of the forces that act on both masses.
            % 
            % % Evaluate the forces acting on $m_t_1$. Define the equations of motion of 
            % the first bob by balancing the horizontal and vertical force components. Specify 
            % these two equations as symbolic equations |eqx_1| and |eqy_1|.

            % eqx_t_1 = M_1*ax_t_1(t) == -T_1*sin(theta_t_1(t)) + T_2*sin(theta_t_2(t))
            % eqy_t_1 = M_1*ay_t_1(t) == +T_1*cos(theta_t_1(t)) - T_2*cos(theta_t_2(t)) - M_1*g
            eqx_t_1 = M_1*ax_t_1(t) == -mo_1*diff(diff(u_1))*cos(theta_t_1(t))-T_1*sin(theta_t_1(t)) + T_2*sin(theta_t_2(t));
            eqy_t_1 = M_1*ay_t_1(t) == -mo_1*diff(diff(u_1))*sin(theta_t_1(t))+T_1*cos(theta_t_1(t)) - T_2*cos(theta_t_2(t)) - M_1*g;
            %% 
            % % Evaluate the forces acting on $m_t_2$. Define the equations of motion of 
            % the second bob by balancing the horizontal and vertical force components. Specify 
            % these two equations as symbolic equations |eqx_2| and |eqy_2|.

            eqx_t_2 = M_2*ax_t_2(t) == -mo_2*diff(diff(u_2))*cos(theta_t_2(t))-T_2*sin(theta_t_2(t));
            eqy_t_2 = M_2*ay_t_2(t) == -mo_2*diff(diff(u_2))*sin(theta_t_2(t))+T_2*cos(theta_t_2(t)) - M_2*g;
            %% Step 3: Evaluate Forces and Reduce System Equations
            % Four equations of motion describe the kinematics of the double pendulum. Evaluate 
            % the forces acting on the rods and reduce the set of four equations to two equations.
            % 
            % The equations of motion have four unknowns: $\theta_1$, $\theta_2$, $T_1$, 
            % and $T_2$. Evaluate the two unknowns $T_1$ and $T_2$ from |eqx_1| and |eqy_1|. 
            % Use <docid:symbolic_ug#buezrr6 solve> function to find $T_1$ and $T_2$.

            Tension = solve([eqx_t_1 eqy_t_1],[T_1 T_2]);
            %% 
            % Substitute the solutions for $T_1$ and $T_2$ into |eqx_2| and |eqy_2|.

            eqRed_1 = subs(eqx_t_2,[T_1 T_2],[Tension.T_1 Tension.T_2]);
            eqRed_2 = subs(eqy_t_2,[T_1 T_2],[Tension.T_1 Tension.T_2]);

%             [Vsym,Ssym] = odeToVectorField(eqRed_1,eqRed_2); % Not necessary
        end
        
        function [sols,xs] = solveEOMs(obj,eqRed_1,eqRed_2)
%             eqn_1 = subs(eqRed_1); % without any further arguments, simply substitutes that which has been defined
%             eqn_2 = subs(eqRed_2);
            eqn_1 = eqRed_1; % without any further arguments, simply substitutes that which has been defined
            eqn_2 = eqRed_2;
            %% 
            % The two equations are nonlinear second-order differential equations. To solve 
            % these equations, convert them to first-order differential equations by using 
            % the <docid:symbolic_ug#bvlo61d odeToVectorField> function.

            [V,S] = odeToVectorField(eqn_1,eqn_2); % V appears to contain the actual symbolic equations (RHS of the 1st order eqn, i.e. the first derivative of the state), and S contains state variable names (LHS)

            %% 
            % The elements of the vector |V| represent the first-order differential equations 
            % that are equal to the time derivative of the elements of |S|. The elements of  
            % |S| are the state variables $\theta_2$, $d\theta_2 /\mathrm{dt}$, $\theta_1$, 
            % and $d\theta_1 /\mathrm{dt}$. The state variables describe the angular displacements 
            % and velocities of the double pendulum.

            %% 
            % Next, convert the first order-differential equations to a MATLAB function 
            % with the handle |M|.

            M = matlabFunction(V,'vars',{'t','Y'}); % Var Y is replacing state vector here
            %% 
            % Define the initial conditions of the state variables as |[pi/4 0 pi/6 0]|. 
            % Use the <docid:matlab_ref#bu00_4l ode45> function to solve for the state variables. 
            % The solutions are a function of time within the interval |[0 10]|.

            sols = ode45(M,obj.trange,obj.X0);
            xs = obj.As*sin(2*pi*obj.fs*sols.x); % History of inputs
        end
        
        function plotsols(obj,sols,input, ANIM_SOL, SAVE_GIF)
            if nargin < 4
                ANIM_SOL = false;
            end
            if nargin < 5
                SAVE_GIF = false;
            end
            
            % Plot the solutions of the state variables.
            xs = input; % Assumes input = xs
            figure()
            plot(sols.x,sols.y) % sols is a struct with x being the independent var and y being the state variables
            hold on
            plot(sols.x,xs)
            legend('\theta_t_2','d\theta_t_2/dt','\theta_t_1','d\theta_t_1/dt','x_susPoint')
            title('Solutions of State Variables')
            xlabel('Time (s)')
            ylabel('Solutions')
            
            if ANIM_SOL
                % Create the animation of the oscillating double pendulum.
                % 
                % First, create four functions that use <docid:matlab_ref#bu7iw_j deval> to 
                % evaluate the coordinates of both pendulums from the solutions |sols|.
                As = obj.As;
                fs = obj.fs;
                Ao1 = obj.Ao1;
                fo1 = obj.fo1;
                Ao2 = obj.Ao2;
                fo2 = obj.fo2;
                g = obj.g;
    %             f0_1 = obj.f0_1;
    %             f0_2 = obj.f0_2;
                L_1 = obj.L_1;
                L_2 = obj.L_2;
                mo_1 = obj.mo_1;
                mt_1 = obj.mt_1;
                mo_2 = obj.mo_2;
                mt_2 = obj.mt_2;
                M_1 = obj.M_1;
                M_2 = obj.M_2;
                th_t_1_sol = @(t) deval(sols,t,3);
                th_t_2_sol = @(t) deval(sols,t,1);

                x_0 = @(t) As*sin(2*pi*fs*t); % No need for deval since independent var
                % x_0 = @(t) x_0(t);
                y_0 = @(t) 0;
                x_t_1 = @(t) x_0(t) + L_1*sin(th_t_1_sol(t));
                y_t_1 = @(t) y_0(t) -L_1*cos(th_t_1_sol(t));
                % x_2 = @(t) L_1*sin(deval(sols,t,3))+L_2*sin(deval(sols,t,1));
                x_t_2 = @(t) x_t_1(t)+L_2*sin(th_t_2_sol(t));
                % y_2 = @(t) -L_1*cos(deval(sols,t,3))-L_2*cos(deval(sols,t,1));
                y_t_2 = @(t) y_t_1(t)-L_2*cos(th_t_2_sol(t));

                u_1 = @(t) Ao1*(sin(2*pi*fo1*t)+1);
                u_2 = @(t) Ao2*(sin(2*pi*fo2*t)+1);

                x_o_1 = @(t) x_t_1(t) + u_1(t).*cos(th_t_1_sol(t));
                y_o_1 = @(t) y_t_1(t) + u_1(t).*sin(th_t_1_sol(t));
                x_o_2 = @(t) x_t_2(t) + u_2(t).*cos(th_t_2_sol(t));
                y_o_2 = @(t) y_t_2(t) + u_2(t).*sin(th_t_2_sol(t)); % Uncomment when u isn't constant

                % x_o_1 = @(t) x_t_1(t) + u_1*cos(th_t_1_sol(t));
                % y_o_1 = @(t) y_t_1(t) + u_1*sin(th_t_1_sol(t));
                % x_o_2 = @(t) x_t_2(t) + u_2*cos(th_t_2_sol(t));
                % y_o_2 = @(t) y_t_2(t) + u_2*sin(th_t_2_sol(t));

                %% 
                % Next, create a stop-motion animation object of the first pendulum bob by using 
                % the <docid:symbolic_ug#mw_50c42eff-e177-477e-ac88-e209c15b0819 fanimator> function. 
                % By default, |fanimator| creates an animation object with 10 generated frames 
                % per unit time within the range of |t| from 0 to 10. Plot the coordinates by 
                % using the |plot| function. Set the _x_-axis and _y_-axis to be equal length.

                figure()
                fanimator(@(t) plot(x_0(t),y_0(t),'ks','MarkerSize',10,'MarkerFaceColor','k'));
                hold on;
                fanimator(@(t) plot(x_t_1(t),y_t_1(t),'rs','MarkerSize',3*mt_1,'MarkerFaceColor','r'));
                axis equal;
                fanimator(@(t) plot(x_o_1(t),y_o_1(t),'ro','MarkerSize',3*mo_1,'MarkerFaceColor','r'));
                %% 
                % Next, add the animation objects of the first rigid rod, the second pendulum 
                % bob, and the second rigid rod.


                fanimator(@(t) plot([x_0(t) x_t_1(t)],[y_0(t) y_t_1(t)],'r-'));
                fanimator(@(t) plot([x_t_1(t) x_o_1(t)],[y_t_1(t) y_o_1(t)],'r-'));
                fanimator(@(t) plot(x_t_2(t),y_t_2(t),'gs','MarkerSize',3*mt_2,'MarkerFaceColor','g'));
                fanimator(@(t) plot(x_o_2(t),y_o_2(t),'go','MarkerSize',3*mo_2,'MarkerFaceColor','g'));
                fanimator(@(t) plot([x_t_1(t) x_t_2(t)],[y_t_1(t) y_t_2(t)],'g-'));
                fanimator(@(t) plot([x_o_2(t) x_t_2(t)],[y_o_2(t) y_t_2(t)],'g-'));
                %% 
                % Add a piece of text to count the elapsed time by using the <docid:matlab_ref#f68-481090 
                % text> function. Use <docid:matlab_ref#btfaj9t-1 num2str> to convert the time 
                % parameter to a string.

                fanimator(@(t) text(0.6,0.6,"Timer: "+num2str(t,2))); %%,'AnimationRange',[0 20]
                hold off;
                playAnimation
                if SAVE_GIF
                    writeAnimation("demos_all_DoublePendulum.gif");
                end
            end
        end
        
        function Xc = theta2disp(obj,thetas_thetaDots,xs)
            % Assumes input is:
            %    thetas_thetaDots=
            %               theta_t_2
            %              Dtheta_t_2
            %               theta_t_1
            %              Dtheta_t_1
            theta1 = thetas_thetaDots(3,:);
            theta2 = thetas_thetaDots(1,:);
            Xc(1,:) = xs;
            Xc(2,:) = Xc(1,:) + obj.L_1*sin(theta1);
            Xc(3,:) = Xc(2,:) + obj.L_2*sin(theta2);
        end
        
        function M = M(obj)
            M_1 = obj.M_1;
             M_2 = obj.M_2;
             
             L_1 = obj.L_1;
             L_2 = obj.L_2;
             
             M = [(M_1 + M_2)*L_1^2, M_2*L_1*L_2;
                  M_2*L_1*L_2,       M_2*L_2^2];
        end
        
        function K = K(obj)
            g = obj.g;
            M_1 = obj.M_1;
             M_2 = obj.M_2;
             
             L_1 = obj.L_1;
             L_2 = obj.L_2;
             K = g*[(M_1 + M_2)*L_1, 0;
                    0,               M_2*L_2];
        end
        
        function T = T(obj) % Vector that transforms acceleration of suspension point input
             M_1 = obj.M_1;
             M_2 = obj.M_2;
             
             L_1 = obj.L_1;
             L_2 = obj.L_2;            
            
             Ts = [(M_1 + M_2)*L_1; M_2*L_2];
             To1 = Ts;
             To2 = [M_1*L_1; M_2*L_2];
             
             T = [Ts,To1,To2];
        end
        
        function A = Ajaclin1(obj,X,U,fs) % First order 
             % Note this is not really a jacobian, as it's pulled from the
             % small angle linearized model
             
             dt = obj.dt; % Define**
             ws = 2*pi*fs;
             ncont = 3; % Should == size(x)(1)
             ndisc = 2*ncont;
             M = obj.M();
             K = obj.K();
             T = obj.T();
            L_1 = obj.L_1;
            L_2 = obj.L_2;
             
             L = [ 1     0     0;
                  -1/L_1 1/L_1 0;
                   0    -1/L_2 1/L_2];
             Mtil = [1, zeros(1,2);
                     T, M]*L;
                     
             Ktil = [ws^2, zeros(1,2);
                     zeros(2,1), K]*L;
             
%              if obj.order == 3 % Specify others, this may change with order of system
               Ginv = inv(eye(ncont)+(dt^2)*inv(Mtil)*Ktil);
               A = [2*Ginv,     -Ginv;
                     eye(ncont), zeros(ncont)];             
%              end
        end
        
        function Xtp1 = flin(obj,Xt,Ut,fs)
            % Temporarily using linearized model here
            A = obj.Ajaclin(Xt,Ut,fs);
            
            Xtp1 = A*Xt; % + B*Ut
        end
    end
end

%     %% Step 6. Show real vs measured motion of optic
%         % Define discrete sampling params
%         % NOTE: Math prefers L>Fs
%         t0 = 0;
%         t1 = 10;
%         Fs = 1000;
%         Ts =(1/Fs);
% 
%         ts = t0:Ts:t1;
%         L = numel(ts);
% 
%         x_o_1_real = x_o_1(ts);
%         y_o_1_real = y_o_1(ts);
%         x_o_2_real = x_o_2(ts);
%         y_o_2_real = y_o_2(ts);
% 
%         % Generate FFT of motions
%         f = Fs*((-L/2):(L/2))/L;
%         f = f(1:end-1); % Find better way to do this*
%         xh_o_1_real = fft(x_o_1_real);
%         yh_o_1_real = fft(y_o_1_real);
%         xh_o_2_real = fft(x_o_2_real);
%         yh_o_2_real = fft(y_o_2_real);
% 
%         % Convert to half-spectra
%         % i.e. Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
%         fhalf = Fs*(0:(L/2))/L;
%         xh_o_1_real_A = abs(xh_o_1_real/L);
%         xh_o_1_real_AH = xh_o_1_real_A(1:L/2+1);
%         xh_o_1_real_AH(2:end-1) = 2*xh_o_1_real_AH(2:end-1);
% 
%         % Introduce OSEM as a measurement device
%         OSEM_f3dB = 1; % Hz... guesstimate
%         OSEM_ENBW = 1.57*OSEM_f3dB; % Approximating as a 1st order lowpass
%         OSEM_gain = HAM_SUS_OSEM_noise('L', fhalf); % Small angle, assume only dealing with length
%         % function motion_spec = HAM_SUS_OSEM_noise(dof, freq)
% 
%         % Take amplitude measurement
%         xh_o_1_OSEM_AH = xh_o_1_real_AH.*OSEM_gain;
% 
%         % Plot
%         figure(fgind)
%         fgind = fgind + 1;
%         plot(fhalf,xh_o_1_real_AH,'k-.')
%         hold on
%         plot(fhalf,xh_o_1_OSEM_AH,'b-.')
%         end
