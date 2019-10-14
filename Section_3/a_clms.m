function [y_hat, e, coeffs_h, coeffs_g] = a_clms(y, w, mu, order)
    samp_len = length(w);
    coeffs_h = complex(zeros(order+1,samp_len));
    coeffs_g = complex(zeros(order+1,samp_len));
    y_hat = complex(zeros(1,samp_len));
    e = complex(zeros(samp_len,1));
    for i = 1:samp_len
      x_current = input_generator(w,order,i);
      y_hat(i)= coeffs_h(:,i)'*x_current + coeffs_g(:,i)'*conj(x_current);     
      e(i) = y(i) - y_hat(i);  
      coeffs_h(:,i+1) = coeffs_h(:,i) + mu*conj(e(i))*x_current;
      coeffs_g(:,i+1) = coeffs_g(:,i) + mu*conj(e(i))*conj(x_current);
    end
end