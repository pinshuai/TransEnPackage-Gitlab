function  Tcrit = transen_crit( M, lag, alpha, numiter, nbins)
%TRANSEN_CRIT Find the critical value of the transfer entropy statistic
%that needs to be exceeded for statistical signficance
%   M = matrix of unshifted variables, e.g., [X Y] for calculating the X>Y
%   transfer entropy. lag = time lag. alpha = significance level. numiter =
%   number of Monte Carlo shufflings to perform. nbins = number of bins to
%   use to discretize the probability distributions.
%   Modified by Laurel Larsen 10/28/16.

%Set default parameters if not specified in the function call.
if nargin< 5
    if nargin < 4
        numiter = 1000;
        if nargin < 2
            alpha = 0.05;
        end
    end
end
        

Tss = NaN(1,numiter); %Initialize the distribution of randomly shuffled transfer entropies.

for ii = 1:numiter
%     Mss = shuffle(M);
    [Tss(ii), ~] = transenshuffle(M, lag, nbins); %Calculate transfer entropy for each Monte Carlo shuffling.
end

Tss = sort(Tss);
Tcrit = Tss(round((1-alpha)*numiter)); %Critical transfer entropy

end

