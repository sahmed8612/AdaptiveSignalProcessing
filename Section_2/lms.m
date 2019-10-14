function [y_hat, e, coefficents] = lms(x, z, mu, order)
    samp_len = length(x);
    coefficents = zeros(order, samp_len-1);
    y_hat = zeros(samp_len, 1);
    e = zeros(samp_len, 1);
    for i = order+1:samp_len
        x_hat = x(i:-1:i-order+1);
        y_hat(i) = (coefficents(:, i-order).') * x_hat;
        e(i) = z(i) - y_hat(i);
        coefficents(:, i-order+1) = coefficents(:, i-order) + mu * e(i) * x_hat;
    end
    coefficents=coefficents';
end

