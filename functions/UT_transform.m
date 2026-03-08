function [SigmaP,W] = UT_transform(NoSig,X,Qxx,gamma,n)
%   function to derive the deterministic Sigma points needed in the
%   Unscented Kalman Filter (UKF) and their corresponding weights
%
%   INPUT:
%       NoSig: number of the Sigma points that should be generated
%       X: state vector
%       Qxx: VCM of the states
%       gamma: scale factor in UKF
%       n: number of the unknown parameters (size of the state vector)
%
%   OUTPUT:
%       SigmaP: the generated Sigma points with [n x NoSig] dimension
%       W: weight of the generated Sigma points with [NoSig x 1] dimension

    % Sigma Points
    SigmaP = zeros(n,NoSig);

    SigmaP(:,1) = X;

    Q_cloud = chol((n + gamma) * Qxx);

    for i = 1:1:n
        SigmaP(:, 1+i) = X + Q_cloud(i,:)' ;
        SigmaP(:, 1+i+n) = X - Q_cloud(i,:)' ;
    end

    % Weights
    W = zeros(NoSig,1);

    W(1,1) = gamma/(2*(n + gamma));

    for i = 2:NoSig
        W(i,1) = (1/(2*(n + gamma)));
    end

end
