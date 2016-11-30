function M_lagged = LagData( M_unlagged, shift )
%LagData.m Shift two time-series.
%   M_unlagged is a matrix [X Y..n], where X and Y are column vectors of the
%   variables to be compared. shift is a row vector that says how much each
%   variable in M_unlagged is to be shifted by.

newlength = size(M_unlagged,1)-max(shift)+min(shift);
M_lagged = NaN(newlength, size(M_unlagged,2));
for ii = 1:size(M_lagged,2)
    M_lagged(:,ii) = M_unlagged(shift(ii)-min(shift)+1:size(M_unlagged,1)-max(shift)+shift(ii), ii);
end

end

