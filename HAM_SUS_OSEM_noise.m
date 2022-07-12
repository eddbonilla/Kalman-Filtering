function motion_spec = HAM_SUS_OSEM_noise(dof, freq)
% motion_spec = HAM_SUS_OSEM_noise(dof, freq)
% dof is either ''L'' or ''P'' or ''Y'' one ''one''
% hand picked, based on the data from alog 11686
% use the damping off, HSTS for LF and RT
% slightly better than original spec at 1-10 Hz,
% matches data.
% the ''one'' is noise of one osem, based on the T1 sensor
% BTL Sept 2020

% make sure the output is a column vector.
[rows, ~] = size(freq);
if rows == 1
    freq = freq';
end

one_osem_freq  = [ 0.01,  0.1,   0.5,     1,       10,    20,    100]';
one_osem_data  = [ 1e-8,  1e-9, 2.5e-10, 1.5e-10, 5e-11, 4e-11, 4e-11]';
% original_data = [ NaN,   1e-8, 3e-10, 1.5e-10, 6e-11, 6e-11];
% figure; loglog(one_osem_freq, one_osem_data, ...
%     one_osem_freq, original_data); title('noise model one osem')

logreqnoise   = interp1(log10(one_osem_freq),log10(one_osem_data),log10(freq));
one_osem_spec = 10.^logreqnoise;
% add the 0.01 and update the 0.1 based on the top sensors; 
% this matches the original spec at 0.5 & 1, 
% but is 4e-11 instead of 6e-11 at 1 & 10 Hz. 
% BTL thinks this is a better match.
if strncmpi(dof,'one',1)
    motion_spec  = one_osem_spec;

elseif strncmpi(dof,'L',1)
    motion_spec  = one_osem_spec/sqrt(2);
    
elseif strncmpi(dof,'P',1)
    motion_spec  = one_osem_spec*(7e-9/3e-10);  %23.3, scale the 0.5 Hz

elseif strncmpi(dof,'Y',1)
    motion_spec  = one_osem_spec*(5.2e-10/6e-11);  %8.7, scale at 100 Hz
else
    error('call with either ''L'' or ''P''')
end
end
