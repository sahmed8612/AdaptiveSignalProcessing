function w=noise_generator(N,noise_power)
    w=randn(N, 1);
    w=w-mean(w);
    w=w*sqrt(noise_power)./sqrt(var(w));
end