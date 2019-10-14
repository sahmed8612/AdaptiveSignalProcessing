function [y_hat, e, coeffecients_w]=dft_clms(y, mu, K, gamma)
    samp_len=length(y);
    coeffecients_w=zeros(K,samp_len);
    k=0:K-1;
    y_hat=zeros(1,samp_len-1);
    e=zeros(1,samp_len-1);
    for i=1:samp_len-1
        x=(1/K)*exp(1j*2*i*pi*k/K).';
        y_hat(i)=coeffecients_w(:,i)'*x;
        e(i)=y(i)-y_hat(i);
        coeffecients_w(:,i+1)=(1-gamma*mu)*coeffecients_w(:,i) + mu*conj(e(i))*x;  
    end
end