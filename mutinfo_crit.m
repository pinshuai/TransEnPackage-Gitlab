function MIcrit = mutinfo_crit( M, nbins, alpha, numiter)
%MUTINFO_CRIT Find the critical value of the mutual information
%statistic that needs to be exceeded for statistical significance.
%   M is the matrix where columns are the individual variables and rows are
%   the binned values in time. NaNs are acceptable. nbins are the number of bins. 
%   alpha is the significance level. 
%   numiter is an optional input that dictates the number of Monte
%   Carlo simulations, or shufflings, to implement

if nargin < 4
    numiter = 1000;
    if nargin <3
    alpha = 0.05;
    end
end

    %First, initialize variables
    MIss = NaN(1,numiter);
    for ii = 1:numiter
        Mss = shuffle(M);
        MIss(ii) = mutinfo(Mss,nbins); %Distribution of random mutual information statistics
    end
MIss = sort(MIss); %Sorted vector of mutual information from low to high
MIcrit = MIss(round((1-alpha)*numiter));
end

