% Function: Generalised LMS Algorithm
% Inputs:
%  1. X - input matrix
%  2. y - true output
%  3. alpha - step size
%  4. phi - non-linearity 
% Ouputs:
%  1. y_hat - estimated output
%  2. w - weights of model
%  3. e - absolute error over time


function [y_hat, w, e] = GLMS(X,y,alpha,phi,w_init)
    if (nargin < 4)
        phi = @(x) x;
        w = zeros(size(X,2),1);
    elseif (nargin < 5)
        w = zeros(size(X,2),1);
    elseif(nargin == 5)
        w = w_init;
    end
    
    %Initalise
    e = zeros(size(X,1),1);
    y_hat = zeros(size(X,1), 1);
    
    for i = 1:size(X,1)
        x = X(i,:)';
        y_true = y(i);
        y_hat(i) = phi(w'*x);
        e(i) = y_true - y_hat(i);
        w = w + alpha*e(i)*x;
    end
end