function mod = aaCreateModel(data, C, B, A)
%AACREATEMODEL Create a model from the input polynomials
%   Initial guess is always zeros, the free variables are those indexes
%   which have one

modelinit = idpoly([1 zeros(1,length(A) -1)], [1 zeros(1,length(B)-1)], [1 zeros(1,length(C) -1)]);

modelinit.Structure.a.Free = A;
modelinit.Structure.b.Free = B;
modelinit.Structure.c.Free = C;

mod = pem(data, modelinit);

present(mod)

end

