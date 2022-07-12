classdef sens
    properties
        sensor
        % Options:
%          'ADE_1mm', 'ADE_p25mm', 'Kaman_IPSmeas',     % Displacement Sensors
%          'L4C', 'GS13meas', 'GS13spec',               % 1 Hz Geophones
%          'T240meas', 'T240spec', 'Guralp_40T_Spec',   % Low-freq Seismometers 
%          'PCB_356M132_spec', 'Wilcoxon_797L_spec'     % Accelerometers
%          'OSEM'                                       % Suspension meas
        DOF = 'L'; % Pertinent for sensor = 'OSEM'
    end
    
    methods
        function obj = sens(type)
            obj.sensor = type; % Add some validation step here
        end
        function Yct = measureAtFrequency(obj,fs,Xct,Ut) % Phase out when moving into adaptive filtering**
            % Currently informs filter of disturbance frequency**
            % currently assumes input only at sus point**
            
            % The structure of this function MUST mirror g()*
            if strcmp(obj.sensor,'OSEM')
                p = 2;
                [n,T] = size(Xct);
                
                R = obj.informedNoise(fs,p);
                V = sqrtm(R)*randn(p,T); % Transform standard gaussian into desired measurement noise
                % OSEM length measures (relative) x coordinate of optics
                
                Yct = obj.g(Xct,Ut) + V;
            elseif strcmp(obj.sensor,'GS13meas') % Remove dependence on L_1,L_2 for this sensor**
                p = 1;
                [n,T] = size(Xct);
                R = obj.informedNoise(fs,p);
                V = sqrtm(R)*randn(p,T); % Transform standard gaussian into desired measurement noise
                % OSEM length measures (relative) x coordinate of optics
                
                Yct = obj.g(Xct,Ut) + V;
            else
                print("Sensor not supported")
                Yct = 0;
            end
        end
        
        function R = informedNoise(obj,f,p) % Phase out when moving into adaptive filtering**
            if strcmp(obj.sensor,'OSEM')
                v_f = HAM_SUS_OSEM_noise(obj.DOF, f); % Assumes a single pertinent DoF
            else
                v_f = SEI_sensor_noise(obj.sensor, f); % m/rtHz
            end
            alpha = 1e10; % Scaling factor to units of m^2
            r2 = alpha*v_f; 
            R = r2*eye(p); % Arbitrary structure, maybe change
        end
        
        function Yc_ideal_t = g(obj,Xct,Ut) % U is a dead input
            % Output is noiseless measurement
            % Assumes state vector has the following structure:
            %             Xt =
            %              xt0
            %              xt1
            %              xt2
            %     NOTE: may have more than this number of entries, so long as these are first         
            if strcmp(obj.sensor,'OSEM')
                p = 2;                
                Yc_ideal_t = [Xct(2,:)-Xct(1,:);Xct(3,:)-Xct(1,:)];
            elseif strcmp(obj.sensor,'GS13meas') % Remove dependence on L_1,L_2 for this sensor**
                p = 1;
                Yc_ideal_t = Xct(1,:);
            end
        end
        
        function C = Cjaclin(obj,Xt,Ut)
            % NOTE: Not really a jacobian as the system was pre-linearized
            if strcmp(obj.sensor,'OSEM')
                C = [-1 1 0 0 0 0;
                     -1 0 1 0 0 0];
            elseif strcmp(obj.sensor,'GS13meas') % Remove dependence on L_1,L_2 for this sensor**
                C = [1 0 0 0 0 0];
            end
        end
    end
end
