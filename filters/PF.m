function [X_plus,Qxx_plus] = PF(TS1,TS2,l,X,Qxx,Qww,Phi,NP,ep_first,ep_last)
%   function of the Particle Filter (PF) algorithm
%
%   INPUT:
%       TS1: position of the first total station (TS1)
%       TS2: position of the second total station (TS2)
%       l: matrix of the measured distance and angle values by means of the
%       total stations
%       X: state vector
%       Qxx: VCM of the states
%       Qww: VCM of the process noise
%       Phi: transition matrix
%       NP: number of the particles
%       ep_first: the first epoch number
%       ep_last: the last epoch number
%
%   OUTPUT:
%       X_plus: filtered state vector
%       Qxx_plus: VCM of the filtered states
X_plus = zeros(size(l,1),size(X,1));
Qxx_plus = cell(size(l,1),1);
    
%% generate the particles and store them in a variable called "S" %%
S = mvnrnd(X,Qxx,NP)';
    
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
    
        %% prediction step of the filter %%
        % fill out the "ffunStates" function to derive the predicted sample
        % points (continue with the name of the variable being "S")
        Qww_NP = mvnrnd(zeros(1,4),Qww,NP)';
        S = ffunStates(S,Phi) + Qww_NP ;
        
        %% update step of the filter %%
        % fill out the "hfun" function to estimate the observations based on
        % the most recent derived particles. then, derive the difference
        % between these estimated observations and the real sensor data
        % (measurements). store the differences in a variable called "res".
        res = l(ep,:)' - hfun(TS1,TS2,S);
    
        % estimate the likelihood of the particles based on the derived
        % differences (res) and by using "mvnpdf" MATLAB library and store
        % these likelihoods in a vector called "q"
        q = mvnpdf(res',zeros(1,4),Qll) + 1*10^(-99);
    
        % derive the importance weights by using the derived likelihoods and
        % store the normalized importance weights in a vector called "w_tilde"        
        w_tilde = (q/NP)/sum(q/NP);

        % resample the particles by using the "ResidualResampling" function 
        i_new = ResidualResample(1:NP,w_tilde');
        S = S(:,i_new);
    
      %% derive the state vector "X" along with its VCM "Qxx" %%
        X_plus(ep,:) = mean(S,2);
        Qxx_plus{ep} = cov(S');        
end
end

