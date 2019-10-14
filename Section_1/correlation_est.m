% Using Biased and Unbiased estiamtors of ACF, the correlogram estiamtor function calculate
% the correlogram spectral estimators
function [r,time_diff,Power_spec,Samp_f] = correlation_est(x,bias)
    [r,time_diff]=xcorr(x,bias);
    int = ifftshift(r);
    Power_spec = real(fftshift(fft(int)))./(2*pi);
    Samp_f=time_diff./max(time_diff);
end