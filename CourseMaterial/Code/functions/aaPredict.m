function [y_hat Fk Gk] = aaPredict(y, C, A, k)
%AAPREDICT Summary of this function goes here
%   Detailed explanation goes here
[Fk Gk] = polydivision(C, A, k);

y_hat = filter(Gk,C,y);

end

