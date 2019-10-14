function [y_hat, e, coeffeicents_h] = clms(y, w, mu, order)
    samp_len = length(w);
    coeffeicents_h = complex(zeros(order+1,samp_len));
    y_hat = complex(zeros(1,samp_len));
    e = complex(zeros(samp_len,1)); 
    for i = 1:samp_len
      x_current = input_generator(w,order,i);
      y_hat(i)= coeffeicents_h(:,i)'*x_current ;     
      e(i) = y(i) - y_hat(i);  
      coeffeicents_h(:,i+1) = coeffeicents_h(:,i) + mu*conj(e(i))*x_current;
    end
end