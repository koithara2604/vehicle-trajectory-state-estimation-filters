function [X_plus,Qxx_plus] = EnKF(TS1,TS2,l,X,Qxx,Qww,Phi,N,ep_first,ep_last)
%   function of the Ensemble Kalman Filter (EnKF) algorithm
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
%       N: number of the generated ensembles
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector
%       Qxx_plus: VCM of the filtered states
    X_plus = zeros(size(l,1),size(X,1));
    Qxx_plus = cell(size(l,1),1);
    
    %% generate the ensembles and store them in a variable called "S" %%
    S = mvnrnd(X,Qxx,N)';

for ep = ep_first:ep_last
        %% measurement accuracies
        % store the measurement accuracies in a structure called sigma ->
        % sigma.d1, sigma.a1, sigma.d2, sigma.a2
         sigma.d1 = (40*10^-3) + (20*l(ep,1)*10^-6);
         sigma.a1 = (20*10^-3)*pi/200;
         sigma.d2 = (25*10^-3) + (15*l(ep,3)*10^-6);
         sigma.a2 = (12*10^-3)*pi/200;
                
         %% form the VCM of the measurement accuracies and call it Qll %%
         d1_squared = sigma.d1^2;
         a1_squared = sigma.a1^2;
         d2_squared = sigma.d2^2;
         a2_squared = sigma.a2^2;
         diag_v = [d1_squared, a1_squared, d2_squared, a2_squared];
         Qll = diag(diag_v);
    
        %% prediction step of the filter %%
        % fill out the "ffunStates" function to derive the predicted ensemble
        % points (continue with the name of the variable being "S")
        Qww_N = mvnrnd(zeros(1,4),Qww,N)';
        S = ffunStates(S,Phi) + Qww_N ;
           
        %% derive the mean value of each state from the ensembles and call it XPred %%
        XPred = mean(S,2);
     
        %% derive the error matrix of the ensembles and call it "Ex" %%
        Ex = S - XPred;
    
        % fill out the "hfun" function to estimate the observations based on
        % the most recent derived ensembles. call the predicted observations
        % "SigmaP_PredObv".
        SigmaP_PredObv = hfun(TS1,TS2,S);
    
        % add noise to the derived deterministic predcited observations by
        % taking into account the measurement accuracies (Qll) and using the
        % multivariate normal distribution. continue with the name of the
        % variable being "SigmaP_PredObv".
        noise = mvnrnd(zeros(1,4),Qll,N)';
        SigmaP_PredObv =  SigmaP_PredObv + noise;
    
        % derive the mean of the predicted observations
        PredObv = mean(SigmaP_PredObv,2);
    
        %% derive the error matrix of the observations and call it "El" %%
        El = SigmaP_PredObv - PredObv;
    
        %% derive the VCM of the ensembles and call it Qxx %%
        Qxx = (Ex * Ex')/(N-1);
    
        %% derive the VCM of the observations and call it Qll %%
        Qll = (El * El')/(N-1);
    
        %% derive the cross VCM between the ensembles and observations and call it "Qxl" %%
        Qxl = (Ex * El')/(N-1);
    
        %% derive the Kalman gain and call it "K" %%
        K = Qxl * (Qll^(-1));
    
        %% update the ensembles based on the derived Kalman gain %%
        S = S + (K * (l(ep,:)' - SigmaP_PredObv)); 
    
        %% derive the state vector "X" along with its VCM "Qxx" %%
        X = mean(S,2);
        delta_X = S - X ;
        Qxx = (delta_X * delta_X')/(N-1);
    
        %updating to retun value
        X_plus(ep,:) = X';
        Qxx_plus {ep} = Qxx;        
end
end

