function X = timeseries_to_matrix (timeseries, N)
    n = length(timeseries);
    X = timeseries(hankel(1:n-N+1,n-N+1:n));
    X(end, :) = [];
end

    



