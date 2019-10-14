function [X, y] = Delta_XY(s, Delta, M)
    if ~isvector(s)
        error("input signal must be 1D");
    end
    if ~isscalar(Delta)
        error("delay parameter must be scalar");
    end
    if ~isscalar(M)
        error("lags parameter must be scalar");
    end
    s = reshape(s, [], 1);
    N = size(s, 1);
    s = [zeros(Delta+M-1, 1); s];
    X = zeros(M, N-1);
    y = zeros(1, N-1);
    for n=1:N
        X(:, n) = s(n+(M-1):-1:n);
        y(n) = s(n+(Delta+M-1));
    end
    
end
