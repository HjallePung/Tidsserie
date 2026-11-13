% FUNCTION [C, stdC] = polyProd( A, B, [stdA], [stdB] )
%
% The function computes the polynomial C=conv(A,B), with the corresponding
% std of the coefficients in C (if given the std of the coefficients in the
% polynomials A and B; otherwise stdC returns an empty vector). 
%

% Assuming that the std of the coefficients are independent and Normal-distributed,
% the variance of the C coefficients may be computed as:
%
%   Var(c_k) = \sum_{i+j=k} a_i^2 \sigma_{b_j}^2 + b_j^2 \sigma_{a_i}^2 + 
%                   + \sigma_{a_i}^2 \sigma_{b_j}^2
% 
% where \sigma_{a_k} and \sigma_{b_k} denotes the std of the k:th A and B
% coefficient, respectively. 

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
function [C, stdC] = polyProd( A, B, stdA, stdB )

% Form polynomial product.
C = conv( A, B );

% If given stds of A and B, form the std of C.
if nargin>2
    nA   = length(A)-1;
    nB   = length(B)-1;
    varC = zeros(1, nA+nB+1);
    for k1 = 0:(nA+nB)
        vCk = 0;
        for k2 = 0:nA
            diffk = k1-k2;
            if diffk >= 0 && diffk <= nB
                A_k1 = A(k2+1);
                B_diffk = B(diffk+1);
                vA  = stdA(k2+1)^2;
                vB  = stdB(diffk+1)^2;
                vCk = vCk + vB*A_k1^2 + vA*B_diffk^2 + vA*vB;
            end
        end
        varC(k1+1) = vCk;
    end
    stdC = sqrt(varC);
else
    stdC = [];
end
