function [NP,GRD] = GRID(GRD,SET)
% Code (FORTRAN) obtained from DVD enclosed with 'Convective Heat Transfer' by Cebeci (2002)
% Few changes made (see comments)

if GRD.VGP - 1 <= 0.001
    NP = round(GRD.etaE/GRD.Deta(1) + 1.0001); % added round()
else
    NP = round(log((GRD.etaE/GRD.Deta(1))*(GRD.VGP - 1) + 1)/log(GRD.VGP) + 1.0001); % added round(log) for substitution of alog in FORTRAN
end
GRD.NP = NP; % store original result for case of multiple BL calculations with same grid
if NP <= SET.NPT
    GRD.eta(1) = 0;
    for j = 2:SET.NPT
        GRD.Deta(j,1) = GRD.VGP*GRD.Deta(j-1);        % preallocate
        GRD.A(j,1) = 0.5*GRD.Deta(j-1);               % preallocate
        GRD.eta(j,1) = GRD.eta(j-1) + GRD.Deta(j-1);  % preallocate
    end
else
    fprintf('NP exceeded NPT - Program terminated')
    keyboard
    return % program proceeds in MAIN-file
end