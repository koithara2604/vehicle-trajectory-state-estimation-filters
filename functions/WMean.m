function X_minus = WMean(X,W)
%   function to derive the weighed mean of the predicted Sigma points
%
%   INPUT:
%       X: predcited Sigma points with dimension [n x NoSig] dimension ->
%       n: number of the unknown parameters (size of the state vector) and
%       NoSig: number of the Sigma points in total
%       W: weight of the Sigma points with [NoSig x 1] dimension
%
%   OUTPUT:
%       X_minus: predicted state vector
    
    [m,n] =size(X);
    X_minus = zeros(m,1);
    for i = 1:1:n
        X_minus = X_minus + (W(i,1).*X(:,i));
    end
end