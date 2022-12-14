function [DFDX] = GRADNT(F,X)
% Function GRADNT computes the derivative (gradient) of F(X) with
% respect to X using finite difference techniques
% Requires input to be continuously increasing or decreasing

% Obtained from: NASA TN D-5681
% FORTRAN Program for Calculating compressible laminar and turbulent
% boundary layers in arbitrary pressure gradients
% William D. McNally, Lewis Research Center, Clevelenad, Ohio 44135
% National aeronautics and space administration, May 1970

N1 = length(X) - 1;
for i = 1:N1
    SL(i) = (F(i+1) - F(i))/(X(i+1) - X(i));
    DIST(i) = sqrt((F(i+1) - F(i))^2 + (X(i+1) - X(i))^2);
end
for i = 2:N1
    DFDX(i) = (SL(i)*DIST(i-1) + SL(i-1)*DIST(i))/(DIST(i-1) + DIST(i));
end

DFDX(1) = SL(1) + (SL(2) - SL(3))*DIST(1)/(DIST(1) + DIST(2));
DFDX(N1+1) = SL(N1) + (SL(N1) - SL(N1-1))*DIST(N1)/(DIST(N1) + DIST(N1-1));
