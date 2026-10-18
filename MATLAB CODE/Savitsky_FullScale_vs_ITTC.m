% Savitsky_FullScale_vs_ITTC.m
% Section 5.5: required thrust of the 50 m vessel from the
clear; close all; clc;

outDir = fileparts(mfilename('fullpath'));
if isempty(outDir), outDir = pwd; end
xlsx = fullfile(outDir, 'ITTC1957_50m_Vessel.xlsx');

%% Vessel (model hull of Table 5.1 scaled by SF)
SF  = 50 / 1.6;                       % length scale factor
g   = 9.81;                           % m/s^2
dCf = 0.0004;                         % ATTC roughness allowance, as at model scale

c.name = '50 m vessel, LCG 17.5 m';
c.LCG  = 0.560 * SF;                  % LCG from transom (m)
c.rho  = 1026;                        % sea water density (kg/m^3)
c.nu   = 0.00122 / 1026;              % sea water kinematic viscosity (m^2/s)

hull.b    = 0.400 * SF;               % beam (m)
hull.beta = 18;                       % deadrise (deg)
hull.VCG  = 0.092 * SF;               % VCG above keel (m)
hull.f    = -0.0696 * SF;             % thrust line to CG, negative above CG (m)
hull.eps  = 0;                        % thrust inclination to keel (deg)
hull.Lp   = 1.600 * SF;               % planing length (m)
volS      = 824.9542425380138;        % displaced volume, SF^3 x model volume (m^3)
hull.W    = c.rho * g * volS;         % displacement (N), 8303 kN

Vref = sqrt(g * volS^(1/3));          % denominator of Fn_V (m/s)

%% ITTC 1957 results (ITTC1957_50m_Vessel.xlsx, sheet 'ITTC 1957 50 m vessel')
FnI0 = [0.227265 0.527124 0.828437 1.129107 1.430154 1.728145 2.029979 2.330872 2.629542 2.928002 3.229235 3.530025 3.828457];
knI0 = [4.237405 9.828350 15.446415 21.052476 26.665564 32.221676 37.849458 43.459669 49.028453 54.593313 60.209875 65.818165 71.382518];
RTI0 = [18.643099 148.498753 456.793712 907.026148 1091.644236 1175.668648 1290.863284 1318.157996 1291.988789 1274.930231 1256.624042 1258.341596 1255.523465];
FnI = FnI0;  knI = knI0;  RTI = RTI0;
if isfile(xlsx)
    try
        sh  = 'ITTC 1957 50 m vessel';
        rdF = readmatrix(xlsx, 'Sheet', sh, 'Range', 'D2:D14')';
        rdK = readmatrix(xlsx, 'Sheet', sh, 'Range', 'L2:L14')';
        rdR = readmatrix(xlsx, 'Sheet', sh, 'Range', 'R2:R14')';
        if numel(rdF) == 13 && numel(rdK) == 13 && numel(rdR) == 13 && ...
                all(isfinite([rdF rdK rdR])) && all([rdF rdK rdR] > 0)
            FnI = rdF;  knI = rdK;  RTI = rdR;
        else
            warning('ITTC workbook read incomplete; stored ITTC 1957 values used.');
        end
    catch err
        warning('ITTC workbook not read (%s); stored ITTC 1957 values used.', err.message);
    end
end

%% Savitsky at the 13 test Froude numbers and on a continuous curve
[tauP, DP] = predict(FnI * Vref, c, hull, g, dCf);
MP  = mfactor(FnI, c.LCG, hull);
RS  = DP / 1000;                      % bare Savitsky (kN)
RSM = MP .* DP / 1000;                % Savitsky + M-factor (kN)
use = FnI >= 1.0 & ~isnan(RS);

eS = 100 * (RS  - RTI) ./ RTI;        % difference from ITTC 1957 (%)
eM = 100 * (RSM - RTI) ./ RTI;

FnC = linspace(1.0, 4.0, 121);
[~, DC] = predict(FnC * Vref, c, hull, g, dCf);
RC  = DC / 1000;
RCM = mfactor(FnC, c.LCG, hull) .* RC;

%% Command-window table
fprintf('Required thrust of the 50 m vessel (kN); difference from ITTC 1957 in brackets\n');
fprintf('  Fr     knots   ITTC 1957   Savitsky          + M-factor        P_E ITTC (MW)\n');
for i = 1:numel(FnI)
    if use(i)
        fprintf('  %4.2f   %5.1f   %7.0f     %7.0f (%+5.1f%%)   %7.0f (%+5.1f%%)   %5.1f\n', ...
            FnI(i), knI(i), RTI(i), RS(i), eS(i), RSM(i), eM(i), RTI(i) * knI(i) * 0.514444 / 1000);
    else
        fprintf('  %4.2f   %5.1f   %7.0f        -                  -               %5.1f\n', ...
            FnI(i), knI(i), RTI(i), RTI(i) * knI(i) * 0.514444 / 1000);
    end
end
fprintf('Mean absolute difference (Fr >= 1.0, %d points): Savitsky %.1f%%, + M-factor %.1f%%\n', ...
    sum(use), mean(abs(eS(use))), mean(abs(eM(use))));
[pI, iI] = max(RTI);  [pS, iS] = max(RC);  [pM, iM] = max(RCM);
fprintf('Peak thrust: ITTC %.0f kN at Fr = %.2f; Savitsky %.0f kN at Fr = %.2f; + M-factor %.0f kN at Fr = %.2f\n', ...
    pI, FnI(iI), pS, FnC(iS), pM, FnC(iM));
fprintf('Peak speeds: ITTC %.1f knots; Savitsky %.1f knots; + M-factor %.1f knots\n', ...
    knI(iI), FnC(iS) * Vref / 0.514444, FnC(iM) * Vref / 0.514444);
fprintf('Peak thrust-to-weight ratio (ITTC): %.3f\n', pI * 1000 / hull.W);

%% Model-scale differences at LCG 560 mm (Chapter 5 comparison, same Fr), for panel (b)
e560S = [-10.8 -12.7 -7.8 -9.0 -8.6 -7.3 -7.1 -6.6 -6.7 -6.1];
e560M = [ -8.5  10.4 14.4  6.0  0.4 -2.4 -4.9 -6.1 -7.3 -7.2];

%% Figure: required thrust against Froude number
ext = FnI < 1.19 | FnI > 2.99;        % outside the wetted-area regression range
fig = figure('Color', 'w', 'Position', [100 100 1150 430]);

subplot(1, 2, 1); hold on;
plot(FnI, RTI, 'k-', 'LineWidth', 1.3);
plot(FnI(~ext), RTI(~ext), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
plot(FnI(ext),  RTI(ext),  'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
plot(FnC, RC,  'b--', 'LineWidth', 1.3);
plot(FnC, RCM, 'r-',  'LineWidth', 1.3);
hold off; grid on; box on;
xlim([0 4]); ylim([0 1500]); set(gca, 'XTick', 0:0.5:4);
xlabel('Fr = Fn_\nabla'); ylabel('Required thrust, T = R_{Ts} (kN)'); title('(a)');
legend('ITTC 1957', 'ITTC, regression range', 'ITTC, extrapolated wetted area', ...
       'Savitsky', 'Savitsky + M-factor', 'Location', 'southeast');

subplot(1, 2, 2); hold on;
plot(FnI(use), eS(use), 'b--s', 'LineWidth', 1.3, 'MarkerFaceColor', 'b', 'MarkerSize', 5);
plot(FnI(use), eM(use), 'r-^',  'LineWidth', 1.3, 'MarkerFaceColor', 'r', 'MarkerSize', 5);
plot(FnI(use), e560S, 'b:s', 'MarkerSize', 5);
plot(FnI(use), e560M, 'r:^', 'MarkerSize', 5);
yline(0, 'k-');
hold off; grid on; box on;
xlim([1 4]); set(gca, 'XTick', 1:0.5:4);
xlabel('Fr = Fn_\nabla'); ylabel('Difference from ITTC 1957 / measurement (%)'); title('(b)'); ylim([-15 20]);
legend('Savitsky, 50 m vs ITTC', '+ M-factor, 50 m vs ITTC', ...
       'Savitsky, model vs measured', '+ M-factor, model vs measured', 'Location', 'northeast');

print(fig, '-dpng', '-r200', fullfile(outDir, 'thrust_vs_Fn_50m.png'));

%% Figure: required thrust against speed in knots
knC  = FnC * Vref / 0.514444;         % curve speeds (knots)
fig2 = figure('Color', 'w', 'Position', [100 100 700 430]); hold on;
plot(knI, RTI, 'k-', 'LineWidth', 1.3);
plot(knI(~ext), RTI(~ext), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
plot(knI(ext),  RTI(ext),  'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
plot(knC, RC,  'b--', 'LineWidth', 1.3);
plot(knC, RCM, 'r-',  'LineWidth', 1.3);
hold off; grid on; box on;
xlim([0 75]); ylim([0 1500]); set(gca, 'XTick', 0:10:70);
xlabel('Speed, V_s (knots)'); ylabel('Required thrust, T = R_{Ts} (kN)');
legend('ITTC 1957', 'ITTC, regression range', 'ITTC, extrapolated wetted area', ...
       'Savitsky', 'Savitsky + M-factor', 'Location', 'southeast');
print(fig2, '-dpng', '-r200', fullfile(outDir, 'thrust_vs_knots_50m.png'));
fprintf('\nSaved thrust_vs_Fn_50m.png and thrust_vs_knots_50m.png to %s\n', outDir);

%% Functions (identical to Savitsky_Code_Test_m_Factor_All.m)
function [tau, DW] = predict(V, c, hull, g, dCf)
    % Returns trim (deg) and resistance D (N) here, not D/W
    tau = nan(size(V));
    DW  = nan(size(V));
    for i = 1:numel(V)
        if ~isfinite(V(i)) || V(i) <= 0, continue; end
        res = @(t) planing(t, V(i), c, hull, g, dCf);
        ts  = 1:0.5:15;
        r   = arrayfun(res, ts);
        j   = find(diff(sign(r)) ~= 0, 1);
        if isempty(j), continue; end
        t = fzero(res, ts([j, j + 1]));
        [~, D, lam] = planing(t, V(i), c, hull, g, dCf);
        if lam * hull.b <= hull.Lp && t >= 2
            tau(i) = t;
            DW(i)  = D;
        end
    end
end

function [Mres, D, lam, Vm, Re, Cf, Df] = planing(tau, V, c, hull, g, dCf)
    b    = hull.b;
    W    = hull.W;
    beta = hull.beta;
    Cv   = V / sqrt(g * b);
    CLb  = W / (0.5 * c.rho * V^2 * b^2);
    CL0  = fzero(@(x) x - 0.0065 * beta * x^0.6 - CLb, [CLb, 2 * CLb + 0.1]);
    lam  = fzero(@(x) tau^1.1 * (0.012 * sqrt(x) + 0.0055 * x^2.5 / Cv^2) - CL0, [0.01, 20]);
    CLd  = 0.012 * sqrt(lam) * tau^1.1;
    Vm   = V * sqrt(max(1 - (CLd - 0.0065 * beta * CLd^0.6) / (lam * cosd(tau)), 0.36));
    Re   = Vm * lam * b / c.nu;
    Cf   = schoenherr(Re);
    Df   = (Cf + dCf) * c.rho * Vm^2 * lam * b^2 / (2 * cosd(beta));
    D    = W * tand(tau) + Df / cosd(tau);
    Cp   = 0.75 - 1 / (5.21 * Cv^2 / lam^2 + 2.39);
    lc   = c.LCG - Cp * lam * b;
    a    = hull.VCG - b / 4 * tand(beta);
    Mres = W * ((1 - sind(tau) * sind(tau + hull.eps)) * lc / cosd(tau) - hull.f * sind(tau)) ...
           + Df * (a - hull.f);
end

function M = mfactor(Fn, LCG, hull)
    r = LCG / hull.b;
    M = 0.98 + 2 * r^1.45 * exp(-2 * (Fn - 0.85)) - 3 * r * exp(-3 * (Fn - 0.85));
    M(Fn < 1.0 | LCG / hull.Lp > 0.46) = NaN;
end

function Cf = schoenherr(Re)
    Cf = fzero(@(x) 0.242 / sqrt(x) - log10(Re * x), [1e-4, 0.1]);
end
