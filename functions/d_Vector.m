function d = d_Vector(R,D,X,TS)
%   function to derive the "d" vector related to the constraints based on
%   the Dx = d constraint equation
%
%   INPUT:
%       R: constraint information
%       D: the derived D-Matrix
%       X: state vector
%       TS: position of the total station to which the contraint is related
%
%   OUTPUT:
%       d: d-vector [c x 1] dimension -> c: number of the constraints
         
    d = R - sqrt((X(1)-TS(1))^2 + (X(2)-TS(2))^2) + D*X ;
    
end