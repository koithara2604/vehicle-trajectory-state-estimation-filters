function A = design_matrix(X,TS1,TS2)
%   function to derive the design matrix
%   
%   INPUT:
%       X: state vector
%       TS1: position of the first total station (TS1)
%       TS2: position of the second total station (TS2)
%
%   OUTPUT:
%       A: design matrix with [n x u] dimension -> n: number of the
%       observations and u: number of the unknowns

    h_x = @(X, TS) ((X(1) - TS(1))/sqrt((X(1) - TS(1))^2 + (X(2) - TS(2))^2));
    h_y = @(X, TS) ((X(2) - TS(2))/sqrt((X(1) - TS(1))^2 + (X(2) - TS(2))^2));
    alpha_x = @(X, TS) ((-1*(X(2) - TS(2)))/((X(1) - TS(1))^2 + (X(2) - TS(2))^2));
    alpha_y = @(X, TS) ((X(1) - TS(1))/((X(1) - TS(1))^2 + (X(2) - TS(2))^2));
    
    A = [h_x(X, TS1), h_y(X, TS1), 0, 0;
        alpha_x(X, TS1), alpha_y(X, TS1), 0, 0;
        h_x(X, TS2), h_y(X, TS2), 0, 0;
        alpha_x(X, TS2), alpha_y(X, TS2), 0, 0];
end