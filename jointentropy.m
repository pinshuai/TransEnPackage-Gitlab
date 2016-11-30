%jointentropy.m
% Written by Laurel Larsen on 10/26/16
% M is a two-column matrix that contains the input vectors of data. It may
% contain blanks (NaNs).
% H is the nonnormalized joint entropy
% nvalidpoints is the number of rows (samples) used to calculate the joint
% entropy.
function [H, nvalidpoints]  = jointentropy(M, nbins)
[~, ~, col1cat] = histcounts(M(:,1), nbins); %Which bin the first data column is in
[~,~,col2cat] = histcounts(M(:,2),nbins); %Which bin the second data column is in
col1cat(col2cat==0)=0; %If there is an NaN for any row, assign the other column in that row to the NaN bin too
col2cat(col1cat==0)=0; %See comment above.
jointentcat = (col1cat-1)*nbins+col2cat; %This classifies the joint entropy bin into a number between 1 and nbins^2. 0 is assigned to rows with misisng data.
[N, ~] = histcounts(jointentcat, nbins^2+1, 'BinLimits', [0, nbins^2+1]); %Number of datapoints within each joint entropy bin. Verified.
p = N/sum(N); %Vector of probabilities
H = -sum(p(p>0).*log2(p(p>0)));
nvalidpoints = sum(N);
end