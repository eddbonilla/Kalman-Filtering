function [HAM_asd]=HAM_asd_generator(freq)
load('HAM5_ref_curvesv1.mat')
HAM_asd=interp1(HAM5_ref.freq,HAM5_ref.long,freq).*(2*pi*freq).^2;

end

