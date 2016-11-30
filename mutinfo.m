%MUTINFO.m
% Written by Laurel Larsen on 11/4/16
% M is a two-column matrix that contains the input vectors of data. It may
% contain blanks (NaNs).
% p is the vector of probabilitites associated with each pair of datapoints
% in M. p should be zero for pairs in M with blanks.
% I is actually the percentage of uncertainty reduction in Y (when M is [X Y])
function I = mutinfo(M, nbins)
[~, ~, col1cat] = histcounts(M(:,1), nbins(2)); %Which bin the first data column is in
[~,~,col2cat] = histcounts(M(:,2),nbins(2)); %Which bin the second data column is in
col1cat(col2cat==0)=0; %If there is an NaN for any row, assign the other column in that row to the NaN bin too
col2cat(col1cat==0)=0; %See comment above.
jointentcat = (col1cat-1)*nbins(2)+col2cat; %This classifies the joint entropy bin into a number between 1 and nbins^2. 0 is assigned to rows with misisng data.
[N, ~] = histcounts(jointentcat, nbins(2)^2, 'BinLimits', [1, nbins(2)^2]); %Number of datapoints within each joint entropy bin. Verified.
p = N/sum(N); %Vector of probabilities
[N1, ~, ~] = histcounts(M(:,1), nbins(1)); %Which bin the first data column is in
[N2,~,~] = histcounts(M(:,2),nbins(1)); %Which bin the second data column is in
p1 = N1/sum(N1);
p2 = N2/sum(N2);
Hy = (-sum(p2(p2>0).*log2(p2(p2>0)))); %Shannon entropy of the sink variable. Used to normalize mutual informaiton in the next line. 
I = ((-sum(p1(p1>0).*log2(p1(p1>0)))-sum(p2(p2>0).*log2(p2(p2>0))))+(sum(p(p>0).*log2(p(p>0))))*log2(nbins(1))/log2(nbins(2)))/Hy; %Mutual information, in bits. Joint entropy is scaled to the number of bins in a single dimension.
end