% Sensitivity checks quoted in Section 5.5, using the same functions and
% inputs as Savitsky_Code_Test_m_Factor_All.m.

clear; clc;

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

for k = 1:numel(tc)
    c    = tc(k);
    Vref = sqrt(g * (hull.W / (c.rho * g))^(1/3));
    [tauP, DWP] = predict(c.Fn * Vref, c, hull, g, dCf);
    use = find(c.Fn >= 1.0 & ~isnan(DWP));
    MP  = mfactor(c.Fn, c.LCG, hull);

    eMeasTrim = zeros(size(use)); dTow = eMeasTrim; dFlo = eMeasTrim; dFhi = eMeasTrim;
    for n = 1:numel(use)
        i = use(n);
        V = c.Fn(i) * Vref;
        [~, Dm] = planing(c.tau(i), V, c, hull, g, dCf);
        eMeasTrim(n) = 100 * (Dm / hull.W - c.DW(i)) / c.DW(i);
        res = @(t) planing(t, V, c, setfield(hull, 'eps', -t), g, dCf);
        ts  = 1:0.5:15; r = arrayfun(res, ts); j = find(diff(sign(r)) ~= 0, 1);
        dTow(n) = fzero(res, ts([j, j + 1])) - tauP(i);
        dFlo(n) = predict(V, c, setfield(hull, 'f', hull.f + 0.010), g, dCf) - tauP(i);
        dFhi(n) = predict(V, c, setfield(hull, 'f', hull.f - 0.010), g, dCf) - tauP(i);
    end
    [~, DW0] = predict(c.Fn(use) * Vref, c, hull, g, 0);
    e0S = 100 * (DW0 - c.DW(use)) ./ c.DW(use);
    e0M = 100 * (MP(use) .* DW0 - c.DW(use)) ./ c.DW(use);
    chg = 100 * (DW0 ./ DWP(use) - 1);

    fprintf('\nLCG %s\n', c.name);
    fprintf('  Drag MAE with Savitsky drag evaluated at the measured trim: %.1f%%\n', mean(abs(eMeasTrim)));
    fprintf('  Trim change, horizontal tow instead of parallel to keel: %.3f to %.3f deg\n', min(dTow), max(dTow));
    fprintf('  Trim change, f 10 mm lower / higher: %.3f to %.3f deg / %.3f to %.3f deg\n', min(dFlo), max(dFlo), min(dFhi), max(dFhi));
    fprintf('  Without roughness allowance: bare drag %.1f to %.1f%%; MAE Savitsky %.1f%%, + M-factor %.1f%%\n', ...
        min(chg), max(chg), mean(abs(e0S)), mean(abs(e0M)));
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