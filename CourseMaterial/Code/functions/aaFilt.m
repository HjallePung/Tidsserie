function data = aaFilt(C,A,y)
%AAFILTER filter + cut maxLag first datapoints
data = aaFilter(C, A, y, length(A));
end
