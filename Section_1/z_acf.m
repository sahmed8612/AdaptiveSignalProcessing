function z=z_acf(M,L)
    z=zeros(1,M);
    for i=0:M-1
        z(i+1)=(M-abs(i))/M;
    end
    z=[fliplr(z(2:end)) z zeros(1,L-(2*M-1))];
end