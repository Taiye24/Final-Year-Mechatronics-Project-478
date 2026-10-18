% Savitsky_Code_Test_m_Factor_All.m
% Savitsky (1964) and Savitsky with the Blount and Fox (1976) M-factor:
% predicted drag and trim compared with the towing-tank measurements.

clear; close all; clc;

% Output folder for figures and tables
outDir = fileparts(mfilename('fullpath'));
if isempty(outDir)
    outDir = pwd;
end

hull.b    = 0.400;      % beam (m)
hull.W    = 264.87;     % displacement weight (N)
hull.beta = 18;         % deadrise angle (deg)
hull.VCG  = 0.092;      % VCG above keel (m)
hull.f    = -0.0696;    % thrust line to CG, negative when above CG (m)
hull.eps  = 0;          % thrust inclination to keel (deg)
hull.Lp   = 1.600;      % planing length (m)
g   = 9.81;             % gravitational acceleration (m/s^2)
dCf = 0.0004;           % ATTC standard roughness allowance

% Water properties: ITTC 7.5-02-01-03 (2011) at 15.8 degC (640 mm) and 16.8 degC
% 640 mm, 3 m/s point: Set A corrected in Data_Eval_640.xls (Fn, D/W and tau are the
% means of Sets A-C from Table X.1 of Calculations_Chapters_4_and_5.xlsx)
tc(1) = struct('name', '640 mm', 'LCG', 0.640, 'rho', 998.9783, 'nu', 1.1150e-6, ...
    'Fn',  [0.238107 0.537875 0.839161 1.139841 1.441471 1.738390 2.040486 2.341884 2.643483 2.940551 3.240283 3.540689 3.839838], ...
    'DW',  [0.004009 0.018779 0.062202 0.107947 0.121812 0.130143 0.140295 0.153987 0.163410 0.171417 0.178555 0.187695 0.198168], ...
    'tau', [2.246969 2.324091 2.665680 3.927781 4.687166 4.953500 5.012747 5.332912 5.576606 5.381742 5.081852 4.750394 4.450749]);
tc(2) = struct('name', '560 mm', 'LCG', 0.560, 'rho', 998.8125, 'nu', 1.0867e-6, ...
    'Fn',  [0.227265 0.527124 0.828437 1.129107 1.430154 1.728145 2.029979 2.330872 2.629542 2.928002 3.229235 3.530025 3.828457], ...
    'DW',  [0.002924 0.020507 0.060263 0.117484 0.142915 0.156284 0.173486 0.180134 0.180417 0.182003 0.183856 0.188740 0.193980], ...
    'tau', [3.253203 3.432668 3.957596 5.698560 6.453209 6.784390 7.263353 7.469527 7.055502 6.501091 5.952781 5.472970 5.020689]);
tc(3) = struct('name', '480 mm', 'LCG', 0.480, 'rho', 998.8125, 'nu', 1.0867e-6, ...
    'Fn',  [0.226685 0.526278 0.827478 1.127589 1.425323 1.704229 2.028444 2.330965 2.626603 2.928875 3.228635 3.524176 3.830125], ...
    'DW',  [0.004915 0.022310 0.063445 0.155686 0.183598 0.197312 0.217176 0.211256 0.201704 0.194662 0.192055 0.195439 0.190994], ...
    'tau', [4.401062 4.614271 5.442894 7.691950 8.531165 9.445192 9.876798 9.155522 8.286257 7.509999 6.763408 6.418060 5.684891]);

FnCurve = linspace(1.0, 4.0, 121);

for k = 1:numel(tc)
    c = tc(k);
    Vref = sqrt(g * (hull.W / (c.rho * g))^(1/3));

    [tauC, DWC] = predict(FnCurve * Vref, c, hull, g, dCf);
    MC = mfactor(FnCurve, c.LCG, hull);

    [tauP, DWP] = predict(c.Fn * Vref, c, hull, g, dCf);
    MP = mfactor(c.Fn, c.LCG, hull);
    use = c.Fn >= 1.0 & ~isnan(DWP);

    errS   = 100 * (DWP - c.DW) ./ c.DW;
    errM   = 100 * (MP .* DWP - c.DW) ./ c.DW;
    errTau = 100 * (tauP - c.tau) ./ c.tau;

    fprintf('\nLCG %s\n', c.name);
    fprintf('  Fn     D/W meas  Savitsky  err(%%)  + M-factor  err(%%)  tau meas  tau pred  err(%%)\n');
    for i = find(use)
        fprintf('  %4.2f   %6.4f    %6.4f   %+6.1f    %6.4f    %+6.1f    %5.2f     %5.2f   %+6.1f\n', ...
            c.Fn(i), c.DW(i), DWP(i), errS(i), MP(i) * DWP(i), errM(i), c.tau(i), tauP(i), errTau(i));
    end
    fprintf('  Mean absolute error: Savitsky %.1f%%, + M-factor %.1f%%, trim %.1f%%\n', ...
        mean(abs(errS(use))), mean(abs(errM(use))), mean(abs(errTau(use))));

    cmp = c.Fn >= 1.0;
    FnM = linspace(min(c.Fn(cmp)), max(c.Fn(cmp)), 200);

    figure('Color', 'w', 'Position', [100 100 900 360]);
    subplot(1, 2, 1);
    plot(FnM, pchip(c.Fn(cmp), c.DW(cmp), FnM), 'k-', FnCurve, DWC, 'b--', FnCurve, MC .* DWC, 'r-', 'LineWidth', 1.3);
    hold on; plot(c.Fn(cmp), c.DW(cmp), 'k.', 'MarkerSize', 10); hold off;
    grid on; axis([1 4 0 0.25]); set(gca, 'XTick', 1:0.5:4);
    xlabel('Fn_\nabla'); ylabel('D/\Delta'); title(['(a) Resistance, LCG ' c.name]);
    legend('Measured', 'Savitsky', 'Savitsky + M-factor', 'Location', 'southeast');
    subplot(1, 2, 2);
    plot(FnM, pchip(c.Fn(cmp), c.tau(cmp), FnM), 'k-', FnCurve, tauC, 'b--', 'LineWidth', 1.3);
    hold on; plot(c.Fn(cmp), c.tau(cmp), 'k.', 'MarkerSize', 10); hold off;
    grid on; axis([1 4 0 11]); set(gca, 'XTick', 1:0.5:4);
    xlabel('Fn_\nabla'); ylabel('\tau (deg)'); title(['(b) Trim, LCG ' c.name]);
    legend('Measured', 'Savitsky', 'Location', 'southeast');
    print('-dpng', '-r200', fullfile(outDir, sprintf('savitsky_LCG%d.png', round(1000 * c.LCG))));

    % Error table as PNG
    T = [c.Fn(use)', c.DW(use)', DWP(use)', errS(use)', MP(use)' .* DWP(use)', errM(use)', ...
         c.tau(use)', tauP(use)', errTau(use)'];
    footer = sprintf('Mean absolute error:   Savitsky %.1f%%,   + M-factor %.1f%%,   trim %.1f%%', ...
        mean(abs(errS(use))), mean(abs(errM(use))), mean(abs(errTau(use))));
    saveTable(T, footer, ['Error table, LCG ' c.name], ...
        fullfile(outDir, sprintf('savitsky_LCG%d_table.png', round(1000 * c.LCG))));
end

fprintf('\nFigures and tables saved to: %s\n', outDir);

function saveTable(T, footer, ttl, fname)
    hdr1 = {'Fn_\nabla', 'D/\Delta', 'D/\Delta', 'Error', 'D/\Delta', 'Error', '\tau meas.', '\tau pred.', 'Error'};
    hdr2 = {'', 'measured', 'Savitsky', '(%)', '+ M-factor', '(%)', '(deg)', '(deg)', '(%)'};
    fmts = {'%4.2f', '%6.4f', '%6.4f', '%+.1f', '%6.4f', '%+.1f', '%.2f', '%.2f', '%+.1f'};
    nR = size(T, 1);
    nC = numel(hdr1);
    xc = ((1:nC) - 0.5) / nC;

    fig = figure('Color', 'w', 'Position', [100 100 950 34 * (nR + 5)]);
    ax  = axes('Parent', fig, 'Position', [0 0 1 1], 'XLim', [0 1], 'YLim', [0 nR + 5], ...
               'YDir', 'reverse', 'Visible', 'off');
    hold(ax, 'on');

    text(0.5, 0.7, ttl, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 13);
    for j = 1:nC
        text(xc(j), 1.6, hdr1{j}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
        text(xc(j), 2.2, hdr2{j}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
    end
    plot([0.02 0.98], [2.7 2.7], 'k-', 'LineWidth', 1);
    for i = 1:nR
        for j = 1:nC
            text(xc(j), 3.2 + (i - 1), sprintf(fmts{j}, T(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 10);
        end
    end
    yEnd = 3.2 + nR - 0.5;
    plot([0.02 0.98], [yEnd yEnd], 'k-', 'LineWidth', 1);
    text(0.5, yEnd + 0.8, footer, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);

    print(fig, '-dpng', '-r200', fname);
    close(fig);
end

function [tau, DW] = predict(V, c, hull, g, dCf)
    tau = nan(size(V));
    DW  = nan(size(V));
    for i = 1:numel(V)
        res = @(t) planing(t, V(i), c, hull, g, dCf);
        ts  = 1:0.5:15;
        r   = arrayfun(res, ts);
        j   = find(diff(sign(r)) ~= 0, 1);
        if isempty(j), continue; end
        t = fzero(res, ts([j, j + 1]));
        [~, D, lam] = planing(t, V(i), c, hull, g, dCf);
        if lam * hull.b <= hull.Lp && t >= 2
            tau(i) = t;
            DW(i)  = D / hull.W;
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