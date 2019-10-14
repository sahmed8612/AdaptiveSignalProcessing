function [y, e, coeffs_est] = GNGD_lms(x, z, mu, rho, order)
    samp_len = length(x);
    coeffs_est = zeros(samp_len,order+1);
    y = zeros(samp_len-1, 1);
    e = zeros(samp_len-1, 1);
    reg_fact=ones(samp_len, 1)*1./mu;
    for i = 2:samp_len-1
        x_hat=input_generator(x,order,i);
        y(i) = x_hat'*coeffs_est(i,:)';
        e(i) = z(i) - y(i);
        coeffs_est(i+1,:) = coeffs_est(i,:) + (1 / (reg_fact(i) + x_hat' * x_hat) ) * e(i) * x_hat';        
        x_hat_hat=input_generator(x,order,i-1);
        reg_fact(i+1) = reg_fact(i) + rho * mu * ((e(i) * e(i-1) * x_hat' * x_hat_hat) / (reg_fact(i-1) + x_hat_hat' * x_hat_hat).^2);
    end
    coeffs_est = coeffs_est(:,2:end);
end
