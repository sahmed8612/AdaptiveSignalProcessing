function [x_hat, e, coeffecients] = anc_lms(s, u, mu, order)
    sample_len = length(s);
    coeffecients = zeros(order+1, sample_len);
    x_hat = zeros(1, sample_len);
    e = zeros(sample_len-1, 1);
    for i = 1:sample_len-1
        u_current=input_generator(u,order,i);
        x_hat(i) = coeffecients(:, i)' * u_current;        
        e(i) = s(i) - x_hat(i);
        coeffecients(:, i+1) = coeffecients(:, i) + mu * e(i) * u_current;
    end    
end

