function [x_hat, e, coeffecients] = ale_lms(s, mu, delay, order)
    samp_len = length(s);
    coeffecients = zeros(order+1, samp_len);
    x_hat = zeros(1, samp_len);
    e = zeros(samp_len-1, 1);
    u = [zeros(delay,1); s];
    for i = 1:samp_len-1
        u_current=input_generator(u,order,i);
        x_hat(i) = coeffecients(:, i)' * u_current;        
        e(i) = s(i) - x_hat(i);
        coeffecients(:, i+1) = coeffecients(:, i) + mu * e(i) * u_current;
    end    
end
