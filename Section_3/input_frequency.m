function in_freq = input_frequency(n)
    if n <= 500
        in_freq = 100; 
    elseif n <= 1000
        in_freq = 100 + (n - 500) / 2;
    elseif n <= 1500
        in_freq = 100 + ((n - 1000) / 25)^2;
    else
        in_freq = 0;
    end
end