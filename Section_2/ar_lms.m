function [y_hat, e, coefficents] = ar_lms(x, mu, order)
    sample_len = length(x);
    coefficents = zeros(order, sample_len-1);
    y_hat = zeros(sample_len, 1);
    e = zeros(sample_len, 1);
    for i = order+1:sample_len
        x_hat = x(i-1:-1:i-order);
        y_hat(i) = (coefficents(:, i-order).') * x_hat;
        e(i) = x(i) - y_hat(i);
        coefficents(:, i-order+1) = coefficents(:, i-order) + mu * e(i) * x_hat;
    end
end