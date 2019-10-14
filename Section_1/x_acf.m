function x=x_acf(M,L)
    k=0:1:M-1;
    auto_corr=(M-k)./M;
    x=[auto_corr zeros(1,L-(2*M-1)) fliplr(auto_corr(2:end))];
end