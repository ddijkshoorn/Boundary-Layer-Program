%% PLOTFILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots boundary layer characteristics, after calculation, or from Stored
% Simulation data (after loading stored simulation data)

%% Test cases - NACA0012 Cebeci 2002

% All laminar profiles
% PLT.NSP = [23 28 33 38 43 48 53];

% Plot profiles inside entire transition region
% PLT.NSP = [52 53 54 55 56 57 58 59 60 61 62 63 64];

% All turbulent profiles
% PLT.NSP = [65 70 75 80 85];

% Combination (interesting evolution)
% PLT.NSP = [23 33 43 53 63 73 83 93];

%% Initialize
% close all

% Analyse flow case
if HVR.P2(1) < 1
    disp('Plate Flow')
else
    disp('Stagnation Point Flow')
end

FST = length(SOL); % might differ from NST in case of separation

%% Plot Boundary Layer Integral Properties

% 1) Boundary layer thicknesses
figure
hold on
plot(X(1:FST),BLC.edge)                     % grid
plot(X(1:FST),BLC.delta)                    % velocity
plot(X(1:FST),BLC.delta_ast)                % displacement
plot(X(1:FST),BLC.theta)                    % momentum
plot(X(1:FST),BLC.delta3/2)                 % energy/kinetic; factor 2; plot(X(1:length(BLC.delta4)),BLC.delta)
plot(X(1:FST),BLC.delta4)                   % enthalpy
% plot(X(1:FST),BLC.delta4a)                  % enthalpy (Schlichting)
plot(X(1:FST),BLC.delta4b)                  % total enthalpy loss
plot(X(1:FST),BLC.delta5)                   % entropy
title('Boundary layer thicknesses')
xlabel('$X \, [-]$','Interpreter','latex')
ylabel('\delta [m]')
legend('edge of grid','velocity','displacement','momentum','kinetic','enthalpy','total enthalpy loss','entropy','Location','northwest')

% 2) Form factor
figure
hold on
plot(X(1:FST),BLC.H,'g')
title('Form factor')
xlabel('$X \, [-]$','Interpreter','latex')
ylabel('$H \quad [-]$','Interpreter','latex')

% 3) Skin friction coefficient
figure
hold on
plot(X(1:FST),BLC.Cf)
plot(X(1:FST),BLC.Cf_L)
plot(X(1:FST),BLC.Cf/2)
title('Skin friction')
xlabel('$X \, [-]$','Interpreter','latex')
ylabel('$C_{\mathrm{f}} \, [-]$','Interpreter','latex')
axis([0 1.1*X(length(SOL)) 0 1.1*max(BLC.Cf)])
% legend('Cf aerodynamics','Local Cf','Cf/2 aerodynamics')%,'SU2 simulation')
legend('$C_{\mathrm{f}} = \frac{\tau_{w}}{\frac{1}{2} \rho u_{\infty}^2}$','$C_{f,\mathrm{local}} = \frac{\tau_{w}}{\frac{1}{2} \rho u_{\mathrm{edge}}^2}$','$\frac{C_f}{2}$','Interpreter','latex')

% 4) Loss coefficient (Denton 1993)

% empirical relation (fit for zero PG) (Denton 1993):
% Cd_Denton_Laminar = 0.173*BLC.Re_theta.^-1;
% Cd_Denton_Turbulent = 0.0056*BLC.Re_theta.^(-1/6);
% for comparison with Denton (1993):
ReTheta1 = 10:10:600; % in Denton 1993 Fig. 5 this is 20:10:600 but Cd is different from laminar relation
ReTheta2 = 200:10:5000;
Cd_Denton_Laminar1 = 0.173*ReTheta1.^-1; % in Denton 1993 Fig. 5 a slightly different relation seems to be used
Cd_Denton_Turbulent2 = 0.0056*ReTheta2.^(-1/6);
%

figure
semilogx(ReTheta1,Cd_Denton_Laminar1,'r')
hold on
semilogx(ReTheta2,Cd_Denton_Turbulent2,'k')
semilogx(BLC.Re_theta,BLC.Cd,'g')
% hold on
% plot(BLC.Re_theta,Cd_Denton_Laminar,'r')
% plot(BLC.Re_theta,Cd_Denton_Turbulent,'k')
% plot(BLC.Re_theta,BLC.Cd,'g')
title('Loss coefficient (Denton (1993))','Interpreter','latex')
xlabel('$\mathrm{Re}_{\theta} \quad [-]$','Interpreter','latex')
ylabel('$C_{\mathrm{d}} \, [-]$','Interpreter','latex')
legend('Laminar - empirical','Turbulent - empirical','Simulation result')
% axis([10^1 10^4.1 0 0.01])
axis([10^1 round(1.1*max(BLC.Re_theta)) 0 0.01])

% % NB on loglog-scale:
% figure
% loglog(ReTheta1,Cd_Denton_Laminar1,'r')
% hold on
% loglog(ReTheta2,Cd_Denton_Turbulent2,'k')
% loglog(BLC.Re_theta,BLC.Cd,'g')
% % hold on
% % plot(BLC.Re_theta,Cd_Denton_Laminar,'r')
% % plot(BLC.Re_theta,Cd_Denton_Turbulent,'k')
% % plot(BLC.Re_theta,BLC.Cd,'g')
% title('Loss coefficient (Denton (1993))')
% xlabel('Re_{\theta} [-]')
% ylabel('Cd [-]')
% legend('Laminar - empirical','Turbulent - empirical','Simulation result')
% % axis([10^1 10^4.1 0 0.01])

% 5) Entropy generated (local entropy generation rate)
figure
hold on
plot(X(1:FST),BLC.Sa)
title('Local entropy generation rate','Interpreter','latex')
xlabel('$X [-]$','Interpreter','latex')
ylabel('$\dot{S}_{\mathrm{A}} \quad [J/K/s/m^2]$','Interpreter','latex')

%% Adiabatic/Heat Transfer

% suitable second condition?
Sum_q_wall = 0;
for i = 1:FST
    Sum_q_wall = Sum_q_wall + BLC.q_wall(FST);
end

if OPT.BCEE == 0 || Sum_q_wall < 1 % adiabatic flow
    % plot adiabatic parameters

    % 6) Adiabatic Wall Temperature (when zero heat flux -> adiabatic wall temperature)
    for j = 1:FST
        Taw(j) = FLP{j}.T(1);
    end
    figure
    plot(X(1:FST),Taw);
    title('Adiabatic Wall Temperature')
    xlabel('X [-]')
    ylabel('T [K]')

    % 7) Adiabatic Wall Enthalpy
    for j = 1:FST
        Haw(j) = EDG.HtE(j)*SOL{j}.g(1); % sol.g(1);
    end
    figure
    hold on
    plot(X,FRS.HtI*ones(1,MON.NST));
    plot(X,EDG.HsE);
    plot(X(1:FST),Haw);
    title('Adiabatic Wall Enthalpy')
    xlabel('X [-]')
    ylabel('h [J/kg/K]')
    legend('h_{0,e}','h_e','h_{aw}')

    % 8) Enthalpy recovery factor
    
    % Approximation
    r = NaN(1,FST);
    if OPT.TRME == 4 % fully turbulent
        for j = 1:FST % change 'MON.NST' into 'FST'?
            r(j) = FLP{j}.Pr(end)^(1/3);
        end
    elseif MON.tr > 0 % MON.tr == 3 || MON.tr == 2 || MON.tr == 1
        for j = 1:MON.NTR - 1
            r(j) = FLP{j}.Pr(end)^(1/2);
        end
        for j = MON.NTR:FST
            r(j) = FLP{j}.Pr(end)^(1/3);
        end
    else % MON.tr = 0; OPT.TRME = 0; laminar flow
        for j = 1:FST
            r(j) = FLP{j}.Pr(end)^(1/2);
        end
    end
    
    figure
    hold on
    plot(X(1:FST),BLC.Hrecovery)
    plot(X(1:FST),r)
    plot(X(1:FST),ones(1,FST))
    title('Enthalpy recovery factor (Enthalpy recovery at wall)')
    xlabel('X [-]')
    ylabel('H_{recovered} [-]')
    % legend('Amount of kinetic energy recovered at wall')
    legend('Calculated','Approximation','ideal (incompressible)')

else
    % plot heat transfer parameters

    % 6) Stanton-number
    figure
    plot(X(1:FST),BLC.St_x)
    title('Stanton-number')
    xlabel('X [-]')
    ylabel('St_{x} [-]')

    % 7) Heat flux
    figure
    plot(X(1:FST),BLC.q_wall);
    title('heat flux')
    xlabel('X [-]')
    ylabel('q_{wall} [W/m2/K]')

    % 8) Wall temperature
    for j = 1:FST
        Twall(j) = FLP{j}.T(1);
    end
    figure
    plot(X(1:FST),Twall);
    title('Wall temperature')
    xlabel('X [-]')
    ylabel('T [K]')
end

%% Input
% uE = UE/UI, MaE or psE = PsE/PtI

% 9) Boundary Condition (BC) Boundary Layer Edge (Input)
figure
title('Input BC BL Edge','Interpreter','latex')
hold on
if OPT.BCIE == 1
    if ~isfield(INP,'uE') % ~exist("INP.uE","var")
        plot(X,INP.UE)
        xlabel('$X [-]$','Interpreter','latex')
        ylabel('$U_{\mathrm{e}} \quad [\mathrm{m}/\mathrm{s}]$','Interpreter','latex')
    else
        plot(X,INP.uE)
        xlabel('$X [-]$','Interpreter','latex')
        ylabel('$u_{\mathrm{e}} = \frac{U_{\mathrm{e}}}{U_{\infty}} \quad [-]$','Interpreter','latex')
    end
elseif OPT.BCIE == 2
    plot(X,INP.MaE)
    xlabel('$X [-]$','Interpreter','latex')
    ylabel('$\mathrm{Ma}_{\mathrm{e}} = \frac{U_{\mathrm{e}}}{a_{\mathrm{e}}} \quad [-]$','Interpreter','latex')
elseif OPT.BCIE == 3
    if ~isfield(INP,'psE') % ~exist("INP.psE","var")
        plot(X,INP.PsE)
        xlabel('$X [-]$','Interpreter','latex')
        ylabel('$P_{\mathrm{e}} \quad [\mathrm{Pa}]$','Interpreter','latex')
    else
        plot(X,INP.psE)
        xlabel('$X [-]$','Interpreter','latex')
        ylabel('$p_{\mathrm{e}} = \frac{P_{\mathrm{e}}}{P_{\mathrm{0,\infty}}} \quad [-]$','Interpreter','latex')
    end
else
    fprintf('Boundary Condition Boundary Layer Edge different than expected\n')
end

%% Pressure Gradient parameters

% 10) Pressure Gradient parameters m1, m2 and m3 as function of surface coordinate X [-]
m1 = HVR.P1;
m2 = HVR.P2;
m3 = 2*m1-1-m2;

figure
hold on
plot(X(1:MON.NST),m1)
plot(X(1:MON.NST),m2)
plot(X(1:MON.NST),m3)
title('Pressure Gradient Parameters')
xlabel('X [-]')
ylabel('PG [-]')
legend('m1','m2','m3')

%% Numerics

% 11) Number of iterations per station
figure
plot(MON.ITE)
title('number of iterations per station')
xlabel('station number [-]')
ylabel('number of iterations [-]')

% 12) Number of mesh grid extensions per station
figure
plot(MON.MGE)
title('Number of mesh grid extensions per station')
xlabel('station number [-]')
ylabel('number of mesh grid extensions [-]')

%% Plot Boundary Layer Profiles

if ~isempty(PLT.NSP)

    if length(PLT.NSP) == 1
        NS_title = sprintf(' profile at station %s', num2str(PLT.NSP(1)));
    else
        NS_title = sprintf(' profiles at stations %s', num2str(PLT.NSP(1)));
        for i = 2:length(PLT.NSP)
            NS = PLT.NSP(i);
            NS_title = strcat(NS_title,sprintf(', %s', num2str(NS)));
        end
    end

    % 13) Shear stress
    figure
    fv1 = get(groot,'CurrentFigure');
    hold on
    %     title('Dimensionless shear stress inside BL','Interpreter','latex')
    title(strcat('Shear stress',NS_title),'Interpreter','latex')
    xlabel('$y/\delta$ \, [-]','Interpreter','latex')
    ylabel('$\frac{\tau}{\tau_{\mathrm{wall}}} \quad [-]$','Interpreter','latex')
    axis([0 2 -0.2 1.2]);

    % 14) Enthalpy profile
    figure
    fh1 = get(groot,'CurrentFigure');
    hold on
    %     title('Total enthalpy profile inside BL')
    title(strcat('Total enthalpy',NS_title),'Interpreter','latex')
    xlabel('$y/\delta$ \, [-]','Interpreter','latex')
    ylabel('$\frac{h_0}{h_{0,\mathrm{e}}} \quad [-]$','Interpreter','latex')

    % 15) Temperature profile
    figure
    fT1 = get(groot,'CurrentFigure');
    hold on
    %     title('Static temperature profile inside BL')
    title(strcat('Static temperature',NS_title),'Interpreter','latex')
    xlabel('$y/\delta$ \, [-]','Interpreter','latex')
    ylabel('$T \quad [\mathrm{K}]$','Interpreter','latex')

    % 16) Velocity profiles in normal coordinates
    figure
    fu1 = get(groot,'CurrentFigure');
    hold on

    title(strcat('Velocity',NS_title),'Interpreter','latex')
    xlabel('$y/\delta$ \, [-]','Interpreter','latex')
    ylabel('$u = \frac{U}{U_{\mathrm{e}}} \quad [-]$','Interpreter','latex')

%% Label writer

Turb = 0; % no turbulence
% l = 0; % counter for number of stations to plot, see loop
% ll = length(PLT.NSP); % total amount of stations to plot
ll1 = 0; % number of laminar stations to plot
ll2 = 0; % number of transitional and turbulent stations combined to plot
ll3 = 0; % number of turbulent stations to plot

for l = 1:length(PLT.NSP)

    NS = PLT.NSP(l);

    if TCC.gamma_tr(NS) > 0 % not laminar, meaning fully turbulent or in transition to turbulent flow
        % turbulent and transitional flow (plot together in linear-scale-plots)
        
        Turb = 1;

        if TCC.gamma_tr(NS) == 1
            % fully turbulent flow

            if ll3 == 0 % first fully turbulent station to plot
                NS_title3 = sprintf(' profile at station %s', num2str(NS));
                NS_prev3 = NS;
            else
                if ll3 == 1
                    NS_title3 = sprintf(' profiles at stations %s', num2str(NS_prev3));
                end
                    NS_title3 = strcat(NS_title3,sprintf(', %s', num2str(NS)));
            end

            ll3 = ll3 + 1;

        end

        if ll2 == 0 % first turbulent/transitional station to plot
            NS_title2 = sprintf(' profile at station %s', num2str(NS));
            NS_prev2 = NS;
        else
            if ll2 == 1
                NS_title2 = sprintf(' profiles at stations %s', num2str(NS_prev2));
            end
            NS_title2 = strcat(NS_title2,sprintf(', %s', num2str(NS)));
        end

        ll2 = ll2 + 1;

    else % elseif TCC.gamma_tr(NS) == 0
        % fully laminar flow

        ll1 = ll1 + 1; % not used

    end
end

%% Label writer

    if Turb == 1 %%%%%|| Trans == 1
        % create figures displaying turbulent properties

        % 17) Eddy viscosity
        figure
        fEV1 = get(groot,'CurrentFigure');
        hold on
        title(strcat('Nondimensional Eddy Viscosity',NS_title2),'Interpreter','latex')
        xlabel('$y/\delta$ \, [-]','Interpreter','latex')
        ylabel('$\nu^{+} = \frac{\nu_{\mathrm{T}}}{\nu} \quad [-]$','Interpreter','latex')
        %         axis([0 2 0 140])
        axis([0 2 0 10]) % temporary y_max

        % 18) Turbulent Prandtl-number vs molecular Prandtl-number
        figure
        fPrT1 = get(groot,'CurrentFigure');
        hold on
        title(strcat('Turbulent Prandtl-number vs molecular Prandtl-number',NS_title2),'Interpreter','latex')
        xlabel('$y/\delta$ \, [-]','Interpreter','latex')
        ylabel('$\mathrm{Pr}_{\mathrm{T}} \quad [-]$','Interpreter','latex')
        axis([0 1.5 0.6 1.6])

        if ll3 > 0 % if ll2 > ll3 && ll3 > 0 % if Trans == 1 && ll2 > ll3

            % % % % %         if Turb == 1
            % 19) Velocity profiles in law-of-the-wall coordinates
            figure
            fu2 = get(groot,'CurrentFigure');
            ax_fu2 = gca; % fu2.ax = gca;
            ax_fu2.XScale = 'log';
            ax_fu2.YScale = 'linear';
            axis([10^-1 10^4 0 40])
            hold on
            title(strcat('Velocity law-of-the-wall',NS_title3),'Interpreter','latex')
            xlabel('$y^+ \, [-]$','Interpreter','latex')
            ylabel('$u^+ \quad [-]$','Interpreter','latex')

            % 20) Eddy viscosity in law-of-the-wall coordinates
            figure
            fEV2 = get(groot,'CurrentFigure');
            ax_fEV2 = gca;
            ax_fEV2.XScale = 'log';
            ax_fEV2.YScale = 'linear';
            axis([10^-1 10^4 0 10])
            hold on
            title(strcat('Nondimensional Eddy Viscosity law-of-the-wall',NS_title3),'Interpreter','latex')
            xlabel('$y^+ \, [-]$','Interpreter','latex')
            ylabel('$\nu^{+} = \frac{\nu_{\mathrm{T}}}{\nu} \quad [-]$','Interpreter','latex')

%         if ll3 > 0 % if ll2 > ll3 % if ll2 > ll3 && ll3 > 0 % if Trans == 1 && ll2 > ll3
            %         if Turb == 1
            % 21) Turbulent Prandtl-number vs molecular Prandtl-number in law-of-the-wall coordinates
            figure
            fPrT2 = get(groot,'CurrentFigure');
            ax_fPrT2 = gca;
            ax_fPrT2.XScale = 'log';
            ax_fPrT2.YScale = 'linear';
            axis([10^-1 10^4 0.6 1.6]) % axis([10^0 10^4 0.6 1.2])
            axis([10^-4 10^4 0.6 1.6])
            %             axis([10^-1 10^4 0 10])
            hold on
            title(strcat('Turbulent Prandtl-number law-of-the-wall',NS_title3),'Interpreter','latex')
            xlabel('$y^+ \, [-]$','Interpreter','latex')
            ylabel('$\mathrm{Pr}_{\mathrm{T}} \quad [-]$','Interpreter','latex')
            % % %         end
        end
    end

    max_tau_prev = 0; % setting the right hight for the shear profile plot

    for j = 1:length(PLT.NSP)

        NS = PLT.NSP(j);

        if TCC.gamma_tr(NS) > 0 % if TCC.gamma_tr(NS) == 1
            % calculate properties for turbulence plots % plot turbulent profiles

            % Re-calculate y-coordinate for extended Eddy Viscosity (beyond BL edge)
            CII = 0;
            for k = 2:length(SOL{NS}.c)
                CII = CII + GRD.A(k)*(SOL{NS}.c(k-1) + SOL{NS}.c(k));
            end

            solc = [SOL{NS}.c; ones(length(GRD.eta) - length(SOL{NS}.c),1)];
            flpmu = [FLP{NS}.mu; FLP{NS}.mu(end)*ones(length(GRD.eta) - length(FLP{NS}.mu),1)];

            CI = 0;
            gamma_int = NaN(length(GRD.eta),1);
            for k = 2:length(GRD.eta)
                CI = CI + GRD.A(k)*(solc(k-1) + solc(k));
                gamma_int(k) = 1.0/(1.0 + 5.5*(CI/CII)^6); % function of eta-grid
            end

            %             if exist("EV_max","var")
            %                 EV_prev = EV_max;
            %             end
            EV_max = max(FLP{NS}.EV);
            EV_extended = [FLP{NS}.EV; EV_max*gamma_int(length(FLP{NS}.EV) + 1:end)];
            Pr_extended = [FLP{NS}.Pr; FLP{NS}.Pr(end)*ones(length(GRD.eta) - length(FLP{NS}.Pr),1)];
            PrT_extended = [FLP{NS}.PrT; FLP{NS}.PrT(end)*ones(length(GRD.eta) - length(FLP{NS}.PrT),1)];
            %             PrT_extended = [FLP{NS}.PrT; FLP{NS}.PrT(end)*ones((length(FLP{NS}.EV) + 1:length(EV_extended)))];

            % Re-calculate y_plus for extended EV (beyond BL edge)
            y = sqrt(EDG.rhoE(NS)*EDG.muE(NS)*X(NS)/EDG.UE(NS))*cumtrapz(GRD.eta,solc./EDG.rhoE(NS));
            y_plus = y.*BLC.u_shear(NS)./(flpmu./EDG.rhoE(NS).*solc); % [-]

            % Number of grid-points to plot in extended EV-plot (determined empirically)
            N = round(1.2*length(BLC.y(:,end)));
            if N > length(GRD.eta)
                N = length(GRD.eta);
            end

            if TCC.gamma_tr(NS) == 1
                % Plot figures
                Data_Legend_Name = sprintf('Turbulent @ Station %s',num2str(NS));
                Data_Legend_Name2 = sprintf('Station %s',num2str(NS));

                figure(fv1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS),'DisplayName',Data_Legend_Name)
                max_tau = max(BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS));
                max_tau_vector(j) = max_tau;
                if max_tau > 1.2 && max_tau > max_tau_prev
                    max_tau_prev = max(max_tau_vector);
                    axis([0 2 -0.2 round(1.1*max(BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS)),1)]);
                    % axis([0 2 -0.2 1.2]);
                end

                figure(fh1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),SOL{NS}.g,'DisplayName',Data_Legend_Name)

                figure(fT1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.T/EDG.TsE(NS),'DisplayName',Data_Legend_Name)

                figure(fu1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),SOL{NS}.u,'DisplayName',Data_Legend_Name)

                figure(fu2) % only plotted when profile is turbulent
                semilogx(BLC.y_plus(1:length(SOL{NS}.u),NS),BLC.u_plus(1:length(SOL{NS}.u),NS),'DisplayName',Data_Legend_Name2)
                %                 hold on

                figure(fEV1)
                plot(y(1:N)/BLC.delta(NS),EV_extended(1:N),'DisplayName',Data_Legend_Name2)
                if EV_max > 10
                    axis([0 2 0 round(1.1*EV_max)]);
                end

%                 % Plot grid:
%                 for i = 1:length(GRD.eta) % length(FLP{NS}.EV)
%                     plot([GRD.eta(i),GRD.eta(i)], [0,round(max(FLP{NS}.EV),-2)],'k')
%                 end

                figure(fEV2)
                semilogx(y_plus(1:N),EV_extended(1:N),'DisplayName',Data_Legend_Name2)

%                 % Plot grid:
%                 for i = 1:length(GRD.eta) % length(FLP{NS}.EV)
%                     plot([GRD.eta(i),GRD.eta(i)], [0,round(max(FLP{NS}.EV),-2)],'k')
%                 end

                figure(fPrT1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.PrT,'DisplayName',sprintf('Turbulent Pr-number @ station %s',num2str(NS)))
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.Pr,'DisplayName',sprintf('Molecular Pr-number @ station %s',num2str(NS)))

                figure(fPrT2)
                %             semilogx(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.PrT,'DisplayName',sprintf('Turbulent Pr-number @ station %s',num2str(NS)))
                %             hold on
                %             semilogx(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.Pr,'DisplayName',sprintf('Molecular Pr-number @ station %s',num2str(NS)))

                %                 semilogx(y(1:N)/BLC.delta(NS),PrT_extended(1:N),'DisplayName',Data_Legend_Name2)
                % %                 hold on
                %                 semilogx(y(1:N)/BLC.delta(NS),Pr_extended(1:N),'DisplayName',Data_Legend_Name2)

                semilogx(y(1:N)/BLC.delta(NS),PrT_extended(1:N),'DisplayName',sprintf('Turbulent Pr-number @ station %s',num2str(NS)))
                %                 hold on
                semilogx(y(1:N)/BLC.delta(NS),Pr_extended(1:N),'DisplayName',sprintf('Molecular Pr-number @ station %s',num2str(NS)))


            else
                % plot profile in transition region
                fprintf('Plotted profiles of station %d inside transition region\n',NS)

                Data_Legend_Name = sprintf('Transitional @ Station %s',num2str(NS));

                figure(fv1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS),'DisplayName',Data_Legend_Name)
                max_tau = max(BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS));
                max_tau_vector(j) = max_tau;
                if max_tau > 1.2 && max_tau > max_tau_prev
                    max_tau_prev = max(max_tau_vector);
                    axis([0 2 -0.2 round(1.1*max(BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS)),1)]);
                    % axis([0 2 -0.2 1.2]);
                end

                figure(fh1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),SOL{NS}.g,'DisplayName',Data_Legend_Name)

                figure(fT1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.T/EDG.TsE(NS),'DisplayName',Data_Legend_Name)

                figure(fu1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),SOL{NS}.u,'DisplayName',Data_Legend_Name)

                figure(fEV1)
                plot(y(1:N)/BLC.delta(NS),EV_extended(1:N),'DisplayName',Data_Legend_Name)
                if EV_max > 10
                    axis([0 2 0 round(1.1*EV_max)]);
                end

%                 % Plot grid:
%                 for i = 1:length(GRD.eta) % length(FLP{NS}.EV)
%                     plot([GRD.eta(i),GRD.eta(i)], [0,round(max(FLP{NS}.EV),-2)],'k')
%                 end

                figure(fPrT1)
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.PrT,'DisplayName',sprintf('Turbulent Pr-number @ station %s',num2str(NS)))
                plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.Pr,'DisplayName',sprintf('Molecular Pr-number @ station %s',num2str(NS)))

            end
            %         if Turb == 1
            if exist("EV_max","var")
                EV_maximum(j) = EV_max;
            end

        elseif TCC.gamma_tr(NS) == 0
            % plot laminar profiles

            Data_Legend_Name = sprintf('Laminar @ Station %s',num2str(NS));

            figure(fv1)
            plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS),'DisplayName',Data_Legend_Name)
            max_tau = max(BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS));
            max_tau_vector(j) = max_tau;
            if max_tau > 1.2 && max_tau > max_tau_prev
                max_tau_prev = max(max_tau_vector);
                axis([0 2 -0.2 round(1.1*max(BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS)),1)]);
                % axis([0 2 -0.2 1.2]);
            end

            figure(fh1)
            plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),SOL{NS}.g,'DisplayName',Data_Legend_Name)

            figure(fT1)
            plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),FLP{NS}.T/EDG.TsE(NS),'DisplayName',Data_Legend_Name)

            figure(fu1)
            plot(BLC.y(1:length(SOL{NS}.u),NS)/BLC.delta(NS),SOL{NS}.u,'DisplayName',Data_Legend_Name)

        end
    end

    figure(fv1)
    %     axis([0 2 -0.2 round(1.1*max(BLC.tau(1:length(SOL{NS}.u),NS)/BLC.tau(1,NS)),1)]);
    %     axis([0 2 -0.2 1.2]);
    set(legend, 'Interpreter', 'latex', 'Location', 'northeast');
    legend('show')

    figure(fh1)
    set(legend, 'Interpreter', 'latex', 'Location', 'southeast');
    legend('show')
    ax_fh1 = gca;
    ax_fh1.XLim = [0 2];

    figure(fT1)
    set(legend, 'Interpreter', 'latex', 'Location', 'northeast');
    legend('show')
    ax_fT1 = gca;
    ax_fT1.XLim = [0 2];

    figure(fu1)
    set(legend, 'Interpreter', 'latex', 'Location', 'southeast');
    legend('show')
    ax_fu1 = gca;
    ax_fu1.XLim = [0 2];
    ax_fu1.YLim = [0 1.2];

    if Turb > 0

        figure(fEV1)
        set(legend, 'Interpreter', 'latex', 'Location', 'northeast');
        legend('show')

%         % Plot grid:
%         for i = 1:length(GRD.eta) % length(FLP{NS}.EV)
%             plot([GRD.eta(i),GRD.eta(i)], [0,round(max(FLP{NS}.EV),-2)],'k')
%             % plot([GRD.eta(i),GRD.eta(i)], [0,200],'k') % large number in case of
%             multiple profiles
%         end

        figure(fPrT1)
        set(legend, 'Interpreter', 'latex', 'Location', 'northeast');
        legend('show')

        if ll3 > 0 % if ll2 > ll3 && ll3 > 0  % if Trans == 1 && ll2 > ll3
            figure(fu2)
            set(legend, 'Interpreter', 'latex', 'Location', 'southeast');
            legend('show')
            
            figure(fEV2)
            set(legend, 'Interpreter', 'latex', 'Location', 'northwest');
            legend('show')
            axis([10^0 10^4 0 round(1.1*EV_max)])

%             % Plot grid:
%             for i = 1:length(GRD.eta) % length(FLP{NS}.EV)
%                 plot([GRD.eta(i),GRD.eta(i)], [0,round(max(FLP{NS}.EV),-2)],'k')
%             end

            figure(fPrT2)
            set(legend, 'Interpreter', 'latex', 'Location', 'northeast');
            legend('show')
        end
    end
end

% Clear help variables
clear max_tau max_tau_prev max_tau_vector

%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
