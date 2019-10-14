function [y_hat_sig, e, coeffecients_h] = ar_clms(x, mu, order)
    samp_len = length(x);
    coeffecients_h = zeros(order,samp_len);
    y_hat_sig = zeros(1,samp_len);
    e = zeros(samp_len,1);
    for i = 1:samp_len-1
      x_current = input_generator(x,order-1,i-1);
      y_hat_sig(i)= coeffecients_h(:,i)'*x_current ;     
      e(i) = x(i) - y_hat_sig(i);  
      coeffecients_h(:,i+1) = coeffecients_h(:,i) + mu*conj(e(i))*x_current;
    end
end