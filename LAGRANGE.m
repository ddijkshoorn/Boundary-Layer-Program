function [DFDX] = LAGRANGE(F,X)
% Code (FORTRAN) obtained from DVD enclosed with 'Convective Heat Transfer' by Cebeci (2002)

NXT = length(X);
for i = 2:NXT
    if i==NXT
        A1 = (X(i-1) - X(i-2))*(X(i) - X(i-2));
        A2 = (X(i-1) - X(i-2))*(X(i) - X(i-1));
        A3 = (X(i) - X(i-1))*(X(i) - X(i-2));
        DFDX(i) = (X(i) - X(i-1))/A1*F(i-2) - (X(i) - X(i-2))/A2*F(i-1) + (2.0*X(i) - X(i-2) - X(i-1))/A3*F(i);
    else
        A1 = (X(i) - X(i-1))*(X(i+1) - X(i-1));
        A2 = (X(i) - X(i-1))*(X(i+1) - X(i));
        A3 = (X(i+1) - X(i))*(X(i+1) - X(i-1));
        DFDX(i) = -(X(i+1) - X(i))/A1*F(i-1) + (X(i+1) - 2.0*X(i) + X(i-1))/A2*F(i) + (X(i) - X(i-1))/A3*F(i+1);
    end
end