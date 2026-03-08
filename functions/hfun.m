function h = hfun(TS1,TS2,X)
%   function to estimate the osbervations based on the most recent state
%   vector
%
%   INPUT:
%       TS1: position of the first total station (TS1)
%       TS2: position of the second total station (TS2)
%       X: state vector
%
%   OUTPUT:
%       h: estimated observation vector
    
    d = @(X, TS) (sqrt((X(1,:) - TS(1)).^2 + (X(2,:) - TS(2)).^2));
    alpha = @(X, TS) atan2((X(2,:) - TS(2)),(X(1,:) - TS(1)));
    
    h = [d(X, TS1); alpha(X, TS1); d(X, TS2); alpha(X, TS2)];

end