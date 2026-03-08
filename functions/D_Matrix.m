function D = D_Matrix(X,TS)
%   function to derive the "D" matrix related to the constraints based on
%   the Dx = d constraint equation
%
%   INPUT:
%       X: state vector
%       TS: position of the total station to which the constraint is
%       related
%
%   OUTPUT:
%       D: D-Matrix with [c x u] dimension -> c: number of the constraints
%       and u: number of the unknowns

     h_x = @(X, TS) ((X(1) - TS(1))/sqrt((X(1) - TS(1))^2 + (X(2) - TS(2))^2));
     h_y = @(X, TS) ((X(2) - TS(2))/sqrt((X(1) - TS(1))^2 + (X(2) - TS(2))^2));

     D = [h_x(X, TS), h_y(X, TS), 0, 0];
end