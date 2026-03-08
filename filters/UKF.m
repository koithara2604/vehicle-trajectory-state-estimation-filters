function [X_plus,Qxx_plus] = UKF(TS1,TS2,l,X,Qxx,Qww,Phi,gamma,ep_first,ep_last)
%   function of the Unscented Kalman Filter (UKF) algorithm
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
%       gamma: scale factor in UKF
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector
%       Qxx_plus: VCM of the filtered states
    X_plus = zeros(size(l,1),size(X,1));
    Qxx_plus = cell(size(l,1),1);
    n = length(X);
    NoSig  = (2*n) + 1;
    
    %% calculte the number of the necessary Sigma points and store them in a variable called NoSig %%
    
for ep = ep_first:ep_last
        %% measurement accuracies
        % store the measurement accuracies in a structure called sigma ->
        % sigma.d1, sigma.a1, sigma.d2, sigma.a2
        sigma.d1 = (40*10^-3) + (20*l(ep,1)*10^-6);
        sigma.a1 = (20*10^-3)*pi/200;
        sigma.d2 = (25*10^-3) + (15*l(ep,3)*10^-6);
        sigma.a2 = (12*10^-3)*pi/200;
            
        %% form the VCM of the measurement accuracies and call it Qll %%
        diag_v = [sigma.d1^2, sigma.a1^2, sigma.d2^2, sigma.a2^2];
        Qll = diag(diag_v);
    
        %% Generating the Sigma points of the previous time step %%
        % fill out the "UT_transform" function to derive the Sigma points
        % (SigmaPoints) and their corresponding weights (WSigma)
        [SigmaPoints, WSigma] = UT_transform(NoSig,X,Qxx,gamma,n);
        
        %% prediction step of the filter %%
        % fill out the "ffunStates" function to derive the predicted Sigma points 
        % (predSigmaP)
        predSigmaP = zeros(n,NoSig);
        for i = 1:NoSig
            predSigmaP(:,i) = ffunStates(SigmaPoints(:,i),Phi);
        end
    
        % fill out the "WMean" function to derive the predicted state vector
        % (predX) that is a weighted mean of the predicted Sigma points
        % (predSigmaP)
        predX =  WMean(predSigmaP,WSigma);
        
        delta_X = predSigmaP - predX;
    
        % fill out the "QMean" function to derive the predicted VCM of the
        % states (predQxx)
        predQxx = QMean(delta_X,WSigma) + Qww ;
        
        % derive the new Sigma points based on the predicted states by using
        % the "UT_transform" function and call them newPredSigmaP with new
        % corresponding weights (newWSigmaP)
        [newPredSigmaP,newWSigmaP] = UT_transform(NoSig,predX,predQxx,gamma,n);
        
        % fill out the "hfun" function to estimate the observations (updSigmaP)
        % based on the most recent Sigma point
        updSigmaP =[];
        for i = 1:1:NoSig
            updSigmaP(:,i) = hfun(TS1,TS2,newPredSigmaP(:,i));
        end
    
        % estimate the predicted observation vector based on the estimated
        % observations (updSigmaP) and by using the "WMean" function
        pred_l = WMean(updSigmaP,newWSigmaP);
    
        delta_l = updSigmaP - pred_l;
    
        % estimate the predicted VCM of the observations (predQll) by using the
        % "QMean" function
        predQll = QMean(delta_l,newWSigmaP) + Qll ;
        
        % fill out the "crossQMean" to derive the cross VCM between the state
        % vector and observation vector and call it Qxl
        Qxl = crossQMean(delta_X,delta_l,newWSigmaP) ;
    
        % derive the Kalman gain (K)
        K = Qxl * predQll^(-1);
        
        % derive the updated state vector along with its corresponding VCM
        X = predX + K * (l(ep,:)' - pred_l);
        Qxx = predQxx - (K * predQll * K');
    
        %updating to retun value
        X_plus(ep,:) = X';
        Qxx_plus {ep} = Qxx;    
end
end