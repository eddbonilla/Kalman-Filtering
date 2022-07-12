function [noise_ASD, response_FR, response_ZPK] = SEI_sensor_noise(sensor_name, freq,verbose)
% SEI_sensor_noise  Noise estimates for HAM-ISI, BSC-ISI, and HEPI SEI sensors
%  version 4
%
%  [noise, vargargout] = SEI_sensor_noise('sensor_name', freq, varargin)
%  noise is the ASD of the self noise of the sensor, in m/rtHz
%    this includes the noise of the default LIGO preamp, if any
%    see T0900450 for more details.
%  freq is the frequency vector (Hz)
%  sensor_name is a string defining the amplifier in question
%   so far 
%          'ADE_1mm', 'ADE_p25mm', 'Kaman_IPSmeas',     % Displacement Sensors
%          'L4C', 'GS13meas', 'GS13spec',               % 1 Hz Geophones
%          'T240meas', 'T240spec', 'Guralp_40T_Spec',   % Low-freq Seismometers 
%          'PCB_356M132_spec', 'Wilcoxon_797L_spec'     % Accelerometers
%   have been included as well as 
%          'GS_16bit_ADC' 
%   which returns in Vrms/rtHz, not m/rtHz.
%       (the code is not case sensitive)
%
% This can also be called with an optional second output argument
% which returns the frequency response of the sensitivity of the sensor
% in volts/ meter, e.g.
%
% [noise, response_FR] = SEI_sensor_noise('GS13meas',freq)
% returns the ASD of displacement noise (noise_ASD) for the GeoTech GS13
% 'GS13' as well as the frequency response of the instrument (response_FR)
%
% noise = SEI_sensor_noise('L4C',freq,true)
% returns the ASD of displacement noise (noise_ASD) only, but also spits
% out several comments about the noise when the third boolean argument
% (verbose) is set to true.
%
% 2009-09-24 Authored Brian Lantz, 
% 2010-03-26 ADC noise added from Jay Heefner
% 2011-10-18 BTL update the T240 response
% 2014-02-25 Added factor of two to noise of 1mm CPS ( noise_ASD    = 2*10.^lognoise;)  
% 2014-11-05 Added Guralp  
% 2017-02-03 Added PCB and Wilcoxon Accelerometers, updated help text and
%            sensor input options 
 
  if nargin < 3  %default to not verbose
      verbose = 0;
  end
  
  if verbose
      disp(' ')
      disp('you are using draft version 5. call Brian Lantz for the real one!')
      disp(' ')
  end

if strncmpi(sensor_name,'ADE_1mm',5)
    % see Andy Stein SEI log enty 1311, nov 4 overnight plot
    freq_data   = [.001,  .002,  .01,    .1,    .7,    100];
    noise_data  = [ 5e-8, 9e-9, 2.5e-9, 5e-10, 2e-10, 2e-10];
    
    lognoise     =  interp1(log10(freq_data),log10(noise_data),log10(freq));
    noise_ASD    = 2*10.^lognoise;
    response_FR  = 10/1e-3 * ones(size(freq)); % 10 V per mm
    response_ZPK = tf(10/1e-3);
    
elseif strncmpi(sensor_name,'ADE_p25mm',5)
    freq_data   =   [.001,   .002,    .01,      .1,    .7,    100];
    noise_data  = [1.5e-08, 2.7e-09, 7.5e-10, 1.5e-10, 6e-11, 6e-11];
    
    lognoise    = interp1(log10(freq_data),log10(noise_data),log10(freq));
    noise_ASD   = 10.^lognoise;
    response_FR = 10/0.25e-3 * ones(size(freq)); % 10V/ 0.25 mm
    response_ZPK = tf(10/0.25e-3); 

elseif strncmpi(sensor_name,'Kaman_IPSmeas',5)
    % See E0900426
    freq_data   = [0.001    .01,      .05,   0.1,      0.5,  1,    10,   100 ];
    noise_data  = [1e-7   1.1e-08,   1e-9, 4.5e-10, 3e-10 3e-10, 3e-10, 3e-10];
    
    lognoise    = interp1(log10(freq_data),log10(noise_data),log10(freq));
    noise_ASD   = 10.^lognoise;
    response_FR = 390/1e-3 * ones(size(freq)); % 390 V / mm
    response_ZPK = tf(390/1e-3);
    
elseif strncmpi(sensor_name,'L4C',3)
    freq_data   =   [1E-4 .04,   .52,     .8,       1.4,     4,       10,      100];
    noise_data  = [1E3 1.0e-06, 1.0e-10, 2.3e-11, 7.0e-12, 2.3e-12, 8.5e-13, 8.0e-14];
    
    lognoise    = interp1(log10(freq_data),log10(noise_data),log10(freq));
    noise_ASD   = 10.^lognoise;
    
    response_ZPK = tf(0);  % this is clearly wrong, BTL Oct 18, 2011
    response_FR  = squeeze(freqresp(response_ZPK, 2*pi*freq));
    
elseif strncmpi(sensor_name,'T240spec',5)
    % this is taken directly from the T240_noise_spec_ASD.m code
    data = [...
        0.001	-171.688; ...
        0.003	-179.481; ...
        0.007	-183.636; ...
        0.018	-187.532; ...
        0.044	-189.610; ...
        0.082	-190.130; ...
        0.226	-190.130; ...
        0.530	-189.091; ...
        1.000	-187.273; ...
        2.257	-183.377; ...
        3.634	-180.000; ...
        6.335	-174.286; ...
        9.803	-168.312; ...
        100     -138];        % Last line added by JSK 2012-05-02, Makes 
                              % displacement spectra role of as 1/f^(1/2), 
                              % 'cause that's what it looked like it was 
                              % turning into at 10 Hz. Complete guess, and
                              % it doesn't matter, since we never use the
                              % T240s above 10Hz anyways; it's just to get
                              % noise coupling predictions up to 100 Hz to
                              % not be identically zero
    
    freq_spec = data(:,1);
    w_spec    = 2*pi*freq_spec;
    
    dB_accel_power_spec = data(:,2);
    accel_amp_spec      = 10.^(dB_accel_power_spec/20);
    disp_amp_spec       = accel_amp_spec ./(w_spec.^2);
    lognoise  = interp1(log10(freq_spec), log10(disp_amp_spec), log10(freq));
    noise_ASD = 10.^lognoise;
    
    % original response, changed to resp per manual on Oct 18, 2011, BTL
    %response_FR = 1196 .* 2*pi*freq;  % V/m 
    %disp('warning, the T240 response is only valid between about 20 mHz and 40 Hz')
    % 
    % the new TF comes from pg 10 of the Trillium 240 OBS user guide
    % originally specified in RAD/SEC, with minus signs. 
    %hence the lack of the usual -2*pi*
    spec_k = 2.316e9; % from the manual
    spec_gain = 1196.5;
    
    % this is in units of V/(m/s);
    response_zpk_temp = spec_gain * zpk([0 0 -108 -161], ...
        [-0.01815 + 0.01799i, -0.01815 - 0.01799i, -173, ...
        -196 + 231i, -196 - 231i, -732 + 1415i,  -732 - 1415i], spec_k) ;
    
    % and this is in V/m
    response_ZPK = tf('s') * response_zpk_temp;
    response_FR = squeeze(freqresp(response_ZPK, 2*pi*freq));
    
 
elseif strncmpi(sensor_name,'T240meas',5)
    %eyeballed from the 090310 data sets
    if verbose
        disp('The 0.01 to 0.03 Hz data from T240meas is just 2*T240spec, and not really measured')
        disp('    (the 0.03 Hz point is measured)')
    end
    freq_data   = [1E-4 .01, .03,    .1,    0.3,    1,     3,     10   100];
    noise_data  = [0.01 2*1.4e-7, 2e-8, 1.5e-9, 2e-10, 4e-11, 5e-12, 1.5e-12 1e-13];
     
    lognoise    = interp1(log10(freq_data),log10(noise_data),log10(freq));
    noise_ASD   = 10.^lognoise;

    % original response, changed to resp per manual on Oct 18, 2011, BTL
    %response_FR = 1196 .* 2*pi*freq;  % V/m 
    %disp('warning, the T240 response is only valid between about 20 mHz and 40 Hz')
    % 
    % the new TF comes from pg 10 of the Trillium 240 OBS user guide
    % originally specified in RAD/SEC, with minus signs. 
    %hence the lack of the usual -2*pi*
    spec_k = 2.316e9; % from the manual
    spec_gain = 1196.5;
    
    % this is in units of V/(m/s);
    response_zpk_temp = spec_gain * zpk([0 0 -108 -161], ...
        [-0.01815 + 0.01799i, -0.01815 - 0.01799i, -173, ...
        -196 + 231i, -196 - 231i, -732 + 1415i,  -732 - 1415i], spec_k) ;
    
    % and this is in V/m
    response_ZPK = tf('s') * response_zpk_temp;
    response_FR = squeeze(freqresp(response_ZPK, 2*pi*freq));

elseif strncmpi(sensor_name,'GS13meas',5)
    % taken directly from GS13_noise_measured_March2007 
    % the noise floor of the GS13 measured on the Tech Demo
    %  this is based on the data from the Tech Demo 3/15/2007
    %  using 2 witnesses with good ADCs, the new readouts, and a Q of ~ 5
    %  it is a little bigger than the expected noise based on the
    %  spec sheet for the LT1012 readout.
    %
    %  BTL, Oct 5, 2008
    %
    % see log entry http://ligo.phys.lsu.edu:8080/SEI/1288

    % data from graphclicks
    
    % at low freq, it looks to scale as x10 for every x2 in freq
    % is a power law of -3.32 (expect -3.5 from 1/f and sens scaling)
    % we will use the 3.5 to be conservative
    % 2e-8 * 10^3.5 = 6.3e-5
    % also assume a 1/f falloff at high freq
    full_data = [...
        1E-4   670;...
        0.01   2e-8 * 10^3.5;...
        0.101	2.138e-8;...
        0.201	1.958e-9;...
        0.400	1.892e-10;...
        0.792	1.732e-11;...
        0.994	7.262e-12;...
        1.258	3.887e-12;...
        1.655	2.449e-12;...
        2.588	1.462e-12;...
        3.954	9.725e-13;...
        8.202	4.423e-13;...
        10.710	3.280e-13;...
        24.796	1.137e-13;...
        46.087	5.610e-14;...
        91.208	2.053e-14
        1e3      2e-15;];
    
    data_freq   = full_data(:,1);
    data_noise  = full_data(:,2);
    lognoise   = interp1(log10(data_freq),log10(data_noise),log10(freq));
    noise_ASD  = 10.^lognoise;

    response_ZPK = tf(0);
    response_FR = NaN * ones(size(freq));
    disp(' warning, the freq resp for the GS13 is not yet defined')
    
elseif strncmpi(sensor_name,'GS13calc',6)
    freq_data   =   [.01    .1       .5      .8        1.2     3       10      100];
    noise_data  = [1.9e-5, 6.1e-9, 1.8e-11, 3.0e-12, 1.3e-12, 5.3e-13, 1.4e-13, 1.3e-14];
    
    lognoise    = interp1(log10(freq_data),log10(noise_data),log10(freq));
    noise_ASD   = 10.^lognoise;

    response_ZPK = tf(0);
    response_FR = NaN*ones(size(freq));
    disp(' warning, the freq resp for the GS13 is not yet defined')

elseif strncmpi(sensor_name,'GS_16bit_ADC',3)
    % General Standards 16 bit ADC PCI66-16AI64SSA
    % http://www.generalstandards.com/specs/66_16ai64ssa_c_spec_100710.pdf
    adcGain = 2^16/40; % [ct/V]
     % 0.5 [mHz] to 100 [Hz] data from Keith Riles, pg 44&45 of G1300997.
    % 1000 [Hz] data point is a guess that it stays flat.
     freq_data           = [0.0005,  0.001,   0.01,    0.1,      1,      10,  100, 1000];
    noise_data_ct_rtHz  = [5e-1,   3.9e-1, 9.5e-2, 3.1e-2, 1.05e-2, 6.5e-3, 6e-3, 6e-3]; % [ct/rtHz]
    noise_data_V_rtHz   = noise_data_ct_rtHz / adcGain; % [V/rtHz]
    
    lognoise    = interp1(log10(freq_data),log10( noise_data_V_rtHz),log10(freq));
    noise_ASD   = 10.^lognoise;

    response_ZPK = tf(65536/40);
    response_FR  = squeeze(freqresp(response_ZPK, 2*pi*freq));
    cprintf('-[0.75,0.25 0.25]','          warning, the freq resp for the ADC has units of Volts/rtHz (not m/rtHz)         .\n');
    
elseif strncmpi(sensor_name,'Guralp_40T_Spec',6)
    % Spec taken from data sheet found on paper in PEM drawer at LHO
    % Spec is assumed to be the bottom of the dynamic range, and captured
    % by-eye. I've extrapolated the noise down to 0.01 [Hz] and out to 100
    % [Hz], assuming the slope at either end of the spec does not change.
    period_data         = [0.01 0.02 0.03 0.04 0.07 0.11 0.25    1   10   50   100]; % [s]
    noise_data_accel_dB = [-90  -110 -121 -128 -140 -150 -160 -168 -165 -162  -161]; % [nonsense units]
    
    freq_data  = fliplr(1./period_data);
    bin_width   = (2^(1/6) - 2^(-1/6))*freq_data; % [Hz]. 1/3 octave binwidth, around central frequency of freq_data
    accel_mps2_pkpk = (10.^(fliplr(noise_data_accel_dB)/20)); % [m/s^2, pkpk] "equal to 6 times the rms acceleration for a 1/3rd octave bandwidth"
    accel_mps2_rms_prtHz = accel_mps2_pkpk ./ (6*sqrt(bin_width)); % [(m/s^2)/rtHz]
    noise_data = accel_mps2_rms_prtHz ./ ((2*pi*freq_data).^2); % [m/rtHz]
    
    lognoise   = interp1(log10(freq_data),log10(noise_data),log10(freq));
    noise_ASD  = 10.^(lognoise);
    
    response_ZPK = tf(0);
    response_FR  = squeeze(freqresp(response_ZPK, 2*pi*freq));
    disp(' warning, the freq resp for the Guralp 40T is not yet defined');
    
elseif strncmpi(sensor_name,'PCB_356M132_spec',6)
    % http://www.pcb.com/Linked_Documents/Aerospace/AD_Low_Outgassing_Lowres.pdf
    % PCB Piezotronics model 356M132
    % 500 mV / g = 0.05102 V / (m/s^2)
    % 0.4 Hz to 4 kHz
    
    freq_data = 10.^(0:4);
    accel_ug_prtHz = [18 6 2 0.7 0.2];
    disp_m_prtHz   = 9.8e-6 * accel_ug_prtHz ./ (2*pi*freq_data).^2;
    
    lognoise = interp1(log10(freq_data),log10(disp_m_prtHz),log10(freq));
    noise_ASD = 10.^(lognoise);
    
elseif strncmpi(sensor_name,'Wilcoxon_797L_spec',8)
    % https://buy.wilcoxon.com/amfile/file/download/file_id/85/product_id/429/
    % Meggitt/Wilcoxon model 797L
    % 100 mV/g = 0.010204 V / (m/s^2)
    % 2 Hz to 7 kHz
    
    freq_data = 10.^(1:3);
    accel_ug_prtHz = [8 5 5];
    disp_m_prtHz   = 9.8e-6 * accel_ug_prtHz ./ (2*pi*freq_data).^2;
    
    lognoise = interp1(log10(freq_data),log10(disp_m_prtHz),log10(freq));
    noise_ASD = 10.^(lognoise);
    
else

    
    disp(' error in SEI_sensor_noise  ');
    disp(' defined sensors are: ')
    disp('          Displacement:  ''ADE_1mm'', ''ADE_p25mm'', ''Kaman_IPSmeas'',')
    disp('          1Hz geophone:  ''L4C'', ''GS13meas'', ''GS13spec'',');
    disp(' Low-freq Seismometers:  ''T240meas'', ''T240spec'', ''Guralp_40T_Spec'',')
    disp('        Accelerometers:  ''PCB_356M132_spec'', ''Wilcoxon_797L_spec'',')
    disp('       Data Acquisiton:  ''GS_16bit_ADC''.')
    disp(' still need to define ''STS2'' ')
    noise_ASD   = NaN * ones(size(freq));
    response_FR = NaN * ones(size(freq));
    
end


