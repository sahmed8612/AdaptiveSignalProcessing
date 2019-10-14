function [y, e, coeffs_est] = gass_lms(x, z, initial_mu, rho, order, algo)
    samp_len = length(x);
    coeffs_est = zeros(samp_len,order+1);   
    y = zeros(samp_len-1, 1);    
    e = zeros(samp_len-1, 1);
    mu = ones(samp_len, 1) * initial_mu; 
    phi = zeros(order+1, 1);
   
    for i = 1:samp_len-1
        x_hat=input_generator(x,order,i);
        y(i) = x_hat'*coeffs_est(i,:)';
        e(i) = z(i) - y(i);   
        coeffs_est(i+1,:) = coeffs_est(i,:) + mu(i)*e(i)*x_hat';
        mu(i+1) = mu(i) + rho * e(i) * x_hat' * phi;
        
        switch algo
            case 'ben' 
                phi = (eye(order) - mu(i) * (x_hat)' * x_hat) * phi + e(i) * x_hat;
            case 'ang'
                phi = 0.8 * phi + e(i) * x_hat;
            case 'mat'
                phi = e(i) * x_hat;                 
        end 
    end
    coeffs_est = coeffs_est(:,2:end);
end

