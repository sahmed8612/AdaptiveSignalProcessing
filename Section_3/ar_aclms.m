function [y_hat, e, coeffecients_h, coeffecients_g] = ar_aclms(x, mu, order)
    samp_len = length(x);
    coeffecients_h = zeros(order,samp_len);
    coeffecients_g = zeros(order,samp_len);
    y_hat = zeros(1,samp_len);
    e = zeros(samp_len,1);
    
    for i = 1:samp_len-1
      x_current = input_generator(x,order-1,i-1);
      y_hat(i)= coeffecients_h(:,i)'*x_current + coeffecients_g(:,i)'*conj(x_current);     
      e(i) = x(i) - y_hat(i);  
      coeffecients_h(:,i+1) = coeffecients_h(:,i) + mu*conj(e(i))*x_current;
      coeffecients_g(:,i+1) = coeffecients_g(:,i) + mu*conj(e(i))*conj(x_current);
    end

end