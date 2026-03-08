function [X_plus,Qxx_plus] = IEKF(TS1,TS2,l,X,Qxx,Qww,Phi,trueStates,constraint,noIter,ep_first,ep_last,R)
%   function of the Iterated Extended Kalman Filter (IEKF) algorithm
%
%   INPUT:
%       TS1: position of the first total station (TS1)
%       TS2: positiion of the second total station (TS2)
%       l: matrix of the measured distance and angle values by means of the
%       total stations
%       X: state vector
%       Qxx: VCM of the states
%       Qww: VCM of the process noise
%       Phi: transition matrix
%       trueStates: matrix of the true trajectory values -> [x y constraint]
%       constraint: true -> constraint applies
%                   false -> constraint does not apply
%       noIter: number of the iterations in the update step
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector
%       Qxx_plus: VCM of the filtered states
    X_plus = zeros(size(l,1),size(X,1));
    Qxx_plus = cell(size(l,1),1);
    
    for ep = ep_first:ep_last
       sigma.d1 = (40*10^-3) + (20*l(ep,1)*10^-6);
       sigma.a1 = (20*10^-3)*pi/200;
       sigma.d2 = (25*10^-3) + (15*l(ep,3)*10^-6);
       sigma.a2 = (12*10^-3)*pi/200;
            
       %% form the VCM of the measurement accuracies and call it Qll %%
       diag_v = [sigma.d1^2, sigma.a1^2, sigma.d2^2, sigma.a2^2];
       Qll = diag(diag_v);
        
       %% prediction step of the filter %%
       % fill out the "ffun" function to derive the predicted state vector
       % (predX) and VCM of the states (predQxx)
       [predX,predQxx] = ffun(X,Qxx,Phi,Qww);
        
        %% update step of the filter %%
        X_n = predX;
        for cou = 1:noIter
            % fill out the "design_matrix" function to derive the design matrix
            A = design_matrix(X_n,TS1,TS2);
        
            % calculate the Kalman gain (K)
            K = predQxx * A' * ((A*predQxx*A') + Qll)^(-1);
            
            % fill out the "hfun" function to estimate the observations based on
            % the most recent state vector
            h = hfun(TS1,TS2,X_n);
            
            % estimate the innovation vector (dl)
            dl = (l(ep , :)' - h);
        
            % update the predicted state vector and its corresponding VCM
            X_n = (predX + K*(dl - A * (predX - X_n)));
            Qxx_n = (eye(4) - (K*A)) * predQxx;
        end       
        X = X_n;
    
        % derive the VCM of the estimations and call it "Qxx"
        Qxx = Qxx_n;
        
        % apply the constraints if necessary
        if(constraint && trueStates(ep,3))
            D = D_Matrix(X,TS1);
            d = d_Vector(R,D,X,TS1);
            W = Qxx^(-1);
            X = X - (W^(-1) * D' * (D * W^(-1) *D')^(-1) * ((D * X) -d));
            Qxx = (Qxx - (Qxx * D' * (D * Qxx *D')^(-1) * D* Qxx));
        end

        %updating to retun value
        X_plus(ep,:) = X';
        Qxx_plus {ep} = Qxx;
    end
end