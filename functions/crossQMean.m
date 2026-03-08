function Q_LP = crossQMean(deltaX,deltaL,W)
%   function to derrive the cross VCM 
%
%   INPUT:
%       deltaX: difference bewteen the Sigma points and the mean of the
%       Sigma points, which represents the state vector. dimention is 
%       [n x NoSig] with n: number of the unknown parameters (size of the 
%       state vector) and NoSig: number of the Sigma points in total
%       deltaL: difference between the estimated observations based on the
%       Sigma points and the mean of the estimated observations by the 
%       Sigma points, which represents the observation vector
%       W: weight of the Sigma points with [NoSig x 1] dimension

    [m,n] = size(deltaX);
    Q_LP = zeros(m,m);
    for i = 1:1:n
         Q_LP =  Q_LP + (W(i,1).*deltaX(:,i).*deltaL(:,i)');
    end
end