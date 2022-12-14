function [sol,flp] = FLDPRS(NS,NP,EDG,sol,flp,FLD,OPT)
% Calculates fluid properties inside boundary layer

ITMAX4 = 10; % no need to change beforehand

for j = 1:NP
    
    %% Calorically Perfect Ideal Gas (compressible, incompressible is optional)
    
    if OPT.GASM == 1
        flp.gamma(j,1) = FLD.gamma; % [-]
        flp.Cp(j,1) = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
        if OPT.COMP > 0 % compressible (standard)
            flp.T(j,1) = (EDG.HtE(NS)*sol.g(j) - 0.5*(EDG.UE(NS)*sol.u(j))^2)/flp.Cp(j);
            sol.c(j,1) = flp.T(j)/EDG.TsE(NS);
            if OPT.CPRN < 1 % constant Pr-number
                if OPT.CCRP > 0 % General and variable Chapman-Rubesin parameter (depending on fluid properties)
                    flp.mu(j,1) = FLD.mu_ref*(flp.T(j)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(flp.T(j) + FLD.SLV); % Sutherland's law - visocity
                    flp.C(j,1) = flp.mu(j)/EDG.muE(NS)/sol.c(j);
                else % Predefined (and) constant Chapman-Rubesin parameter
                    flp.mu(j,1) = FLD.C*sol.c(j)*EDG.muE(NS);
                    flp.C(j,1) = FLD.C; % Chapman-Rubesin parameter
                end
                flp.Pr(j,1) = FLD.Pr;
                flp.k(j,1) = flp.mu(j)*flp.Cp(j)/flp.Pr(j);
            else
                flp.mu(j,1) = FLD.mu_ref*(flp.T(j)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(flp.T(j) + FLD.SLV); % Sutherland's law - visocity
                flp.C(j,1) = flp.mu(j)/EDG.muE(NS)/sol.c(j);
                flp.k(j,1) = FLD.k_ref*(flp.T(j)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(flp.T(j) + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
                flp.Pr(j,1) = flp.mu(j)*flp.Cp(j)/flp.k(j);
            end
        else % incompressible, assuming constant fluid properties, only possible in combination with OPT.GASM=1
            sol.c(j,1) = 1;
            flp.C(j,1) = 1;
            flp.T(j,1) = EDG.TsE(NS);
            flp.mu(j,1) = EDG.muE(NS);
            flp.Pr(j,1) = FLD.Pr;
            flp.k(j,1) = flp.mu(j)*flp.Cp(j)/flp.Pr(j);            
        end
        sol.d(j,1) = flp.C(j)*(EDG.UE(NS)^2)*(1 - 1/flp.Pr(j))/EDG.HtE(NS);
        
        %% Thermally Perfect Ideal Gas (compressible by definition)
    
    elseif OPT.GASM == 2

        h0 = EDG.HtE(NS)*sol.g(j);
        U = EDG.UE(NS)*sol.u(j);
        To = 0;
        if j == 1
            Tn = EDG.TtE(NS);
        else
            Tn = flp.T(j - 1,1);
        end
        IT2 = 0;
        while abs(Tn - To) > 1e-3 % probably accurate enough, since method is effective and efficient
            IT2 = IT2 + 1;
            if IT2 > ITMAX4
                fprintf('IT exceeded ITMAX = %g at station = %g \n',ITMAX4, 1)
                keyboard
                return
            end
            To = Tn;
            h = FLD.Cpcoeff(1)*To + FLD.Cpcoeff(2)/2*To^2 + FLD.Cpcoeff(3)/3*To^3 + FLD.Cpcoeff(4)/4*To^4 + FLD.Cpcoeff(5)/5*To^5;
            Cp = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*To + FLD.Cpcoeff(3)*To^2 + FLD.Cpcoeff(4)*To^3 + FLD.Cpcoeff(5)*To^4;
            F = h0 - h - U^2/2;
            dF = -Cp;
            Tn = To - F/dF;
        end

        flp.T(j,1) = Tn;
        flp.Cp(j,1) = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*flp.T(j) + FLD.Cpcoeff(3)*flp.T(j)^2 + FLD.Cpcoeff(4)*flp.T(j)^3 + FLD.Cpcoeff(5)*flp.T(j)^4; % Andrews (1981) IG relation; range: 100 - 590 K
        flp.mu(j,1) = FLD.mu_ref*(flp.T(j)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(flp.T(j) + FLD.SLV); % [Pas], Sutherland's law - visocity
        flp.k(j,1) = FLD.k_ref*(flp.T(j)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(flp.T(j) + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
        flp.Pr(j,1) = flp.mu(j)*flp.Cp(j)/flp.k(j);
        sol.c(j,1) = flp.T(j)/EDG.TsE(NS);
        flp.C(j,1) = flp.mu(j)/EDG.muE(NS)/sol.c(j);
        sol.d(j,1) = flp.C(j)*(EDG.UE(NS)^2)*(1 - 1/flp.Pr(j))/EDG.HtE(NS);
        
        %% Nonideal Gas
    
    elseif OPT.GASM == 3
        flp.T(j,1) = FLD.FP.Temperature('Ph',EDG.PsE(NS)/1e5,(EDG.HtE(NS)*sol.g(j) - (EDG.UE(NS)*sol.u(j))^2/2)/1e3) + 273.15;
        flp.rho(j,1) = FLD.FP.Density('Ph',EDG.PsE(NS)/1e5,(EDG.HtE(NS)*sol.g(j) - (EDG.UE(NS)*sol.u(j))^2/2)/1e3);
        flp.mu(j,1) = FLD.FP.Viscosity('Ph',EDG.PsE(NS)/1e5,(EDG.HtE(NS)*sol.g(j) - (EDG.UE(NS)*sol.u(j))^2/2)/1e3);
        sol.c(j,1) = EDG.rhoE(NS)/flp.rho(j);
        flp.C(j,1) = flp.mu(j)/EDG.muE(NS)/sol.c(j);
        flp.k(j,1) = FLD.FP.ThermCond('Ph',EDG.PsE(NS)/1e5,(EDG.HtE(NS)*sol.g(j) - (EDG.UE(NS)*sol.u(j))^2/2)/1e3);
        flp.Cp(j,1) = FLD.FP.HeatCapP('Ph',EDG.PsE(NS)/1e5,(EDG.HtE(NS)*sol.g(j) - (EDG.UE(NS)*sol.u(j))^2/2)/1e3)*1e3;
        flp.gamma(j,1) = flp.Cp(j)/(flp.Cp(j) - FLD.Rsg);
        flp.Pr(j,1) = flp.mu(j)*flp.Cp(j)/flp.k(j);
        sol.d(j,1) = flp.C(j)*(EDG.UE(NS)^2)*(1 - 1/flp.Pr(j))/EDG.HtE(NS);
        FLD.Z(j,1) = EDG.PsE(NS)/FLD.Rsg/flp.rho(j)/flp.T(j);
        FLD.GamFun(j,1) = FLD.FP.Gamma('Ph',EDG.PsE(NS)/1e5,(EDG.HtE(NS)*sol.g(j) - (EDG.UE(NS)*sol.u(j))^2/2)/1e3);
    end
end
