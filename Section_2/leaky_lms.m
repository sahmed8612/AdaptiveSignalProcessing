function [y_hat, e, coefficents] = leaky_lms(x, mu, gamma, order)
    samp_len = length(x);
    coefficents = zeros(order, samp_len-1);
    y_hat = zeros(samp_len, 1);
    e = zeros(samp_len, 1);
    for i = order+1:samp_len
        x_hat = x(i-1:-1:i-order);
        y_hat(i) = (coefficents(:, i-order).') * x_hat;        
        e(i) = x(i) - y_hat(i);
        coefficents(:, i-order+1) = (1 - mu * gamma).*coefficents(:, i-order) + mu * e(i) * x_hat;
    end
    
end