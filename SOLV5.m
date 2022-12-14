function [sol, DvW] = SOLV5(NP,S,B,R,GRD,HVR,sol)
% Code (FORTRAN) obtained from DVD enclosed with 'Convective Heat Transfer' by Cebeci (2002)
% Few changes made (see comments)
% NB consider preallocating of variables to increase speed

%% Variables to pass
A = GRD.A;
alpha0 = HVR.alpha0;
alpha1 = HVR.alpha1;
% INPUT/PRECAL
% INP: S (use scaled s or S?), eta, eta_inf, alpha0, alpha1, vgp, deta(1), (inside INP!)
% BC: gW or pW, 
% Edge
% Freestream
% Wall

% Calculation: KBM
% sol: f, u, v, g, p, b, c, d, e
% props
% propsE/propsW? freestream? or BC?
% BLC: delta, delta_ast, theta, etc. (BL characteristics)

%% Initialization
% (consider preallocation)
W = zeros(NP,5);
GG = zeros(NP,15);
AA = zeros(NP,15);

%% Elements of triangle-matrix
AA(1,1) = 1;
AA(1,2) = 0;
AA(1,3) = 0;
AA(1,4) = 0;
AA(1,5) = 0;
AA(1,6) = 0;
AA(1,7) = 1;
AA(1,8) = 0;
AA(1,9) = 0;
AA(1,10) = 0;
AA(1,11) = 0;
AA(1,12) = 0;
AA(1,13) = 0;
AA(1,14) = alpha0;
AA(1,15) = alpha1;

%% Elements of W-vector
W(1,1) = R(1,1);
W(1,2) = R(1,2);
W(1,3) = R(1,3);
W(1,4) = R(1,4);
W(1,5) = R(1,5);

%% Forward sweep
% Definitions
for j=2:NP
    AA1=A(j)*AA(j-1,9)-AA(j-1,10);
    AA2=A(j)*AA(j-1,14)-AA(j-1,15);
    AA3=A(j)*AA(j-1,2)-AA(j-1,3);
    AA4=A(j)*AA(j-1,7)-AA(j-1,8);
    AA5=A(j)*AA(j-1,12)-AA(j-1,13);
    AA6=A(j)*AA(j-1,4)-AA(j-1,5);
    AA7=A(j)*S(j,6)-S(j,2);
    AA8=S(j,8)*A(j);
    AA9=A(j)*B(j,6)-B(j,10);
    AA10=A(j)*B(j,8)-B(j,2);
    % Elements of triangle matrix
    DET=AA(j-1,1)*(AA4*AA2-AA1*AA5)-AA(j-1,6)*(AA3*AA2-AA5*AA6)+AA(j-1,11)*(AA3*AA1-AA4*AA6);
    GG(j,1)=(-(AA4*AA2-AA5*AA1)+A(j)^2*(AA(j-1,6)*AA2-AA(j-1,11)*AA1))/DET;
    GG(j,2)=((AA3*AA2-AA5*AA6)-A(j)^2*(AA(j-1,1)*AA2-AA(j-1,11)*AA6))/DET;
    GG(j,3)=(-(AA3*AA1-AA4*AA6)+A(j)^2*(AA(j-1,1)*AA1-AA(j-1,6)*AA6))/DET;
    GG(j,4)=GG(j,1)*AA(j-1,2)+GG(j,2)*AA(j-1,7)+GG(j,3)*AA(j-1,12)+A(j);
    GG(j,5)=GG(j,1)*AA(j-1,4)+GG(j,2)*AA(j-1,9)+GG(j,3)*AA(j-1,14);
    GG(j,6)=(S(j,4)*(AA2*AA4-AA1*AA5)+AA(j-1,11)*(AA1*AA7-AA4*AA8)+AA(j-1,6)*(AA5*AA8-AA7*AA2))/DET;
    GG(j,7)=(AA(j-1,1)*(AA2*AA7-AA5*AA8)+AA(j-1,11)*(AA3*AA8-AA6*AA7)+S(j,4)*(AA5*AA6-AA2*AA3))/DET;
    GG(j,8)=(AA(j-1,1)*(AA4*AA8-AA1*AA7)+S(j,4)*(AA3*AA1-AA4*AA6)+AA(j-1,6)*(AA7*AA6-AA3*AA8))/DET;
    GG(j,9)=GG(j,6)*AA(j-1,2)+GG(j,7)*AA(j-1,7)+GG(j,8)*AA(j-1,12)-S(j,6);
    GG(j,10)=GG(j,6)*AA(j-1,4)+GG(j,7)*AA(j-1,9)+GG(j,8)*AA(j-1,14)-S(j,8);
    GG(j,11)=(B(j,4)*(AA4*AA2-AA5*AA1)-AA9*(AA(j-1,6)*AA2-AA(j-1,11)*AA1)+AA10*(AA(j-1,6)*AA5-AA(j-1,11)*AA4))/DET;
    GG(j,12)=(-B(j,4)*(AA3*AA2-AA5*AA6)+AA9*(AA(j-1,1)*AA2-AA(j-1,11)*AA6)-AA10*(AA(j-1,1)*AA5-AA(j-1,11)*AA3))/DET;
    GG(j,13)=(B(j,4)*(AA3*AA1-AA4*AA6)-AA9*(AA(j-1,1)*AA1-AA(j-1,6)*AA6)+AA10*(AA(j-1,1)*AA4-AA(j-1,6)*AA3))/DET;
    GG(j,14)=GG(j,11)*AA(j-1,2)+GG(j,12)*AA(j-1,7)+GG(j,13)*AA(j-1,12)-B(j,6);
    GG(j,15)=GG(j,11)*AA(j-1,4)+GG(j,12)*AA(j-1,9)+GG(j,13)*AA(j-1,14)-B(j,8);
    % Elements of triangle matrix
    AA(j,1)=1.0;
    AA(j,2)=-A(j)-GG(j,4);
    AA(j,3)=A(j)*GG(j,4);
    AA(j,4)=-GG(j,5);
    AA(j,5)=A(j)*GG(j,5);
    AA(j,6)=S(j,3);
    AA(j,7)=S(j,5)-GG(j,9);
    AA(j,8)=S(j,1)+A(j)*GG(j,9);
    AA(j,9)=-GG(j,10)+S(j,7);
    AA(j,10)=A(j)*GG(j,10);
    AA(j,11)=B(j,3);
    AA(j,12)=B(j,5)-GG(j,14);
    AA(j,13)=B(j,9)+A(j)*GG(j,14);
    AA(j,14)=B(j,7)-GG(j,15);
    AA(j,15)=B(j,1)+A(j)*GG(j,15);
    % Elements of W-vector
    W(j,1)=R(j,1)-GG(j,1)*W(j-1,1)-GG(j,2)*W(j-1,2)-GG(j,3)*W(j-1,3)-GG(j,4)*W(j-1,4)-GG(j,5)*W(j-1,5);
    W(j,2)=R(j,2)-GG(j,6)*W(j-1,1)-GG(j,7)*W(j-1,2)-GG(j,8)*W(j-1,3)-GG(j,9)*W(j-1,4)-GG(j,10)*W(j-1,5);
    W(j,3)=R(j,3)-GG(j,11)*W(j-1,1)-GG(j,12)*W(j-1,2)-GG(j,13)*W(j-1,3)-GG(j,14)*W(j-1,4)-GG(j,15)*W(j-1,5);
    W(j,4)=R(j,4);                                                             
    W(j,5)=R(j,5);
end

%% Backward sweep
j=NP;
% Definitions
DP = -(AA(j,11)*(AA(j,3)*W(j,2)-W(j,1)*AA(j,8))-AA(j,12)*(AA(j,1)*W(j,2)-W(j,1)*AA(j,6)) + W(j,3)*(AA(j,1)*AA(j,8)-AA(j,3)*AA(j,6)));
DV = -(AA(j,11)*(W(j,1)*AA(j,10)-W(j,2)*AA(j,5))-W(j,3)*(AA(j,1)*AA(j,10)-AA(j,5)*AA(j,6))+AA(j,15)*(AA(j,1)*W(j,2)-W(j,1)*AA(j,6)));
DF = -(W(j,3)*(AA(j,3)*AA(j,10)-AA(j,8)*AA(j,5))-AA(j,13)*(W(j,1)*AA(j,10)-AA(j,5)*W(j,2)) + AA(j,15)*(W(j,1)*AA(j,8)-AA(j,3)*W(j,2)));
D1 = -(AA(j,11)*(AA(j,3)*AA(j,10)-AA(j,8)*AA(j,5))-AA(j,13)*(AA(j,1)*AA(j,10)-AA(j,6)*AA(j,5))+AA(j,15)*(AA(j,1)*AA(j,8)-AA(j,6)*AA(j,3)));
% Elements of delta-vector for j=NP
DELP(j) = DP/D1;
DELV(j) = DV/D1;
DELF(j) = DF/D1;
DELG(j) = 0;
DELU(j) = 0;
while j > 1
j = j-1;
% Definitions
BB1=DELU(j+1)-A(j+1)*DELV(j+1)-W(j,4);                                      
BB2=DELG(j+1)-A(j+1)*DELP(j+1)-W(j,5);
CC1=W(j,1)-AA(j,2)*BB1-AA(j,4)*BB2;
CC2=W(j,2)-AA(j,7)*BB1-AA(j,9)*BB2;
CC3=W(j,3)-AA(j,12)*BB1-AA(j,14)*BB2;
DD1=AA(j,3)-AA(j,2)*A(j+1);
DD2=AA(j,8)-AA(j,7)*A(j+1);
DD3=AA(j,13)-AA(j,12)*A(j+1);
EE1=AA(j,5)-AA(j,4)*A(j+1);
EE2=AA(j,10)-AA(j,9)*A(j+1);
EE3=AA(j,15)-AA(j,14)*A(j+1);
DETT=AA(j,1)*DD2*EE3+AA(j,6)*DD3*EE1+AA(j,11)*DD1*EE2-AA(j,11)*DD2*EE1-AA(j,6)*DD1*EE3-AA(j,1)*DD3*EE2;
% Elements of delta-vector
DELF(j)=(CC1*DD2*EE3+CC2*DD3*EE1+CC3*DD1*EE2-CC3*DD2*EE1-CC2*DD1*EE3-CC1*DD3*EE2)/DETT;
DELV(j)=(AA(j,1)*CC2*EE3+AA(j,6)*CC3*EE1+AA(j,11)*CC1*EE2-AA(j,11)*CC2*EE1-AA(j,6)*CC1*EE3-AA(j,1)*CC3*EE2)/DETT;
DELP(j)=(AA(j,1)*CC3*DD2+AA(j,6)*CC1*DD3+AA(j,11)*CC2*DD1-AA(j,11)*CC1*DD2-AA(j,6)*CC3*DD1-AA(j,1)*CC2*DD3)/DETT;
DELU(j)=BB1-A(j+1)*DELV(j);
DELG(j)= BB2-A(j+1)*DELP(j);
end

% New values of F, U, V, G, P
for j=1:NP
    sol.f(j) = sol.f(j) + DELF(j);
    sol.u(j) = sol.u(j) + DELU(j);
    sol.v(j) = sol.v(j) + DELV(j);
    sol.g(j) = sol.g(j) + DELG(j);
    sol.p(j) = sol.p(j) + DELP(j);
end
sol.u(1) = 0;
DvW = DELV(1); % Added
