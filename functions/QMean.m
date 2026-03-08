function Qxx_minus = QMean(delta,W)
%   function to derive the predicted VCM of the states
%
%   INPUT:
%       delta: difference bewteen the Sigma points and the mean of the
%       Sigma points, which represents the state vector. dimention is 
%       [n x NoSig] with n: number of the unknown parameters (size of the 
%       state vector) and NoSig: number of the Sigma points in total
%       W: weight of the Sigma points with [NoSig x 1] dimension
%
%   OUTPUT:
%       Qxx_minus: predicted VCM of the states

    [m,n] = size(delta);
    Qxx_minus = zeros(m,m);
    for i = 1:1:n
        Qxx_minus = Qxx_minus + (W(i,1).*delta(:,i).*delta(:,i)');
    end
    
end