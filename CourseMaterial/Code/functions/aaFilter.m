function data = aaFilter(C,A,y,maxLag)
%AAFILTER filter + cut maxLag first datapoints
data = filter(C,A,y);
data = data(maxLag+1:end);
end

