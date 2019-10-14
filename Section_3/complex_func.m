function y=complex_func(b,a,w)
    samp_len=length(w);
    y=zeros(1,samp_len);
    w=[0 w];
    for i=1:samp_len
        y(i)=a*w(i+1)+b(1)*w(i)+b(2)*conj(w(i));
    end
end
