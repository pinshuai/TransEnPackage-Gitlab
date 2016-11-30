function Mss = shuffle( M )
%SHUFFLE 
%   Shuffles the entries of the matrix M in time while keeping NaNs (blank
%   data values) NaNs. M is a matrix where the columns are individual
%   variables and the rows are entries in time.

Mss = NaN(size(M)); %Initialize
for n = 1:size(M,2) %Break apart columns to shuffle separately
    notnans = find(~isnan(M(:,n))); %Indices of M that are not NaNs.
    R = rand(size(notnans));
    [~,I] = sort(R, 1);
    Mss(notnans, n) = M(notnans(I), n);
end
end

