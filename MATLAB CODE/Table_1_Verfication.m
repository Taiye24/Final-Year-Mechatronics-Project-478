% Savitsky (1964) Table 1 worked example (p. 89): check of the planing
% Chapter 3. This script is the verification scriptT

clear; clc;
hull = struct('b', 14, 'W', 60000, 'beta', 10, 'VCG', 2.0, 'f', 0.50, 'eps', 4, 'Lp', Inf);
c    = struct('LCG', 29.0, 'rho', 1.94, 'nu', 1.0e-5);   % rho from the C_Lb line (Eq. 3.5), nu from Row 6 (Eq. 3.10)
g    = 32.2;       % gravitational acceleration (ft/s^2)
V    = 67.5;       % 40 knots (ft/s)
dCf  = 0.0004;     % ATTC standard roughness allowance, Delta C_A (Eq. 3.12) (Row 8)
tauP  = [2 3 4];   % trial trim angles (deg) used in the paper's Table 1 (trial angles of Eq. 3.18)
pub   = [3.85 2.60 1.86; 67.0 66.6 66.2; 3.61e8 2.42e8 1.73e8; 0.00174 0.00184 0.00192; ...
         7340 5160 3760; 9434 8304 7948; -149960 336600 648350];
% Row labels, with the equation that produces each quantity
names = {'lambda (3.7)', 'Vm (3.8)', 'Re (3.10)', 'Cf (3.11)', 'Df (3.13)', 'D (3.14)', 'Eq. 3.17 resid.'};
res = zeros(7, 3);
for i = 1:3
    [M, D, lam, Vm, Re, Cf, Df] = planing(tauP(i), V, c, hull, g, dCf);
    res(:, i) = [lam; Vm; Re; Cf; Df; D; M];
end
err = 100 * (res - pub) ./ abs(pub);
fprintf('%-16s', 'Row');
fprintf('   tau = %d deg: code, paper, error ', tauP);
fprintf('\n');
for k = 1:7
    fprintf('%-16s', names{k});
    fprintf('   %10.4g %10.4g %+6.1f%%        ', [res(k, :); pub(k, :); err(k, :)]);
    fprintf('\n');
end
% Equilibrium trim by linear interpolation of the Eq. 3.17 residual between
% the bracketing trial angles (2 and 3 deg), Eq. 3.18:
%   tau_eq = tau_i - R_i*(tau_(i+1) - tau_i)/(R_(i+1) - R_i)
te = 2 + res(7, 1) / (res(7, 1) - res(7, 2));
% Drag at equilibrium: the resistance array (Eq. 3.14) is interpolated linearly
% to tau_eq, as described after Eq. 3.18.
De = res(6, 1) + (te - 2) * (res(6, 2) - res(6, 1));
fprintf('\nEquilibrium by linear interpolation between 2 and 3 deg (Eq. 3.18, as in the paper):\n');
fprintf('  trim %.3f deg (paper 2.3 deg, %+.1f%%), drag %.0f lb (paper 9095 lb, %+.1f%%)\n', ...
    te, 100 * (te - 2.3) / 2.3, De, 100 * (De - 9095) / 9095);
% Exact equilibrium: root of the Eq. 3.17 residual, no interpolation
tx = fzero(@(t) planing(t, V, c, hull, g, dCf), [2 3]);
[~, Dx] = planing(tx, V, c, hull, g, dCf);
fprintf('Exact equilibrium (root of Eq. 3.17): trim %.3f deg, drag %.0f lb\n', tx, Dx);
% Paper Row 10 recomputed with Eq. 3.13: Df = (Cf + dCA)*rho*lambda*Vm^2*b^2/(2*cos(beta)),
% where 0.00214 = Cf + dCA (Row 8) at tau = 2 deg
row10 = 0.00214 * 1.94 * 67.0^2 * 3.85 * 14^2 / (2 * cosd(10));
fprintf('\nPaper Row 10 at 2 deg recomputed from its own Rows 4, 5 and 9 with rho = 1.94 (Eq. 3.13): %.0f lb (printed 7340 lb)\n', row10);

function [Mres, D, lam, Vm, Re, Cf, Df] = planing(tau, V, c, hull, g, dCf)
    % Savitsky (1964) planing solution at one trial trim angle tau (deg).
    % Mres is the pitching-moment residual of Eq. 3.17 (Savitsky's Eq. 35).
    % Note: hull.W is the displacement weight (Delta in the document).
    b    = hull.b;
    W    = hull.W;
    beta = hull.beta;
    Cv   = V / sqrt(g * b);                                              % Eq. 3.4  speed coefficient
    CLb  = W / (0.5 * c.rho * V^2 * b^2);                                % Eq. 3.5  deadrise surface load coefficient
    CL0  = fzero(@(x) x - 0.0065 * beta * x^0.6 - CLb, [CLb, 2 * CLb + 0.1]);   % Eq. 3.6  zero-deadrise lift coefficient
    lam  = fzero(@(x) tau^1.1 * (0.012 * sqrt(x) + 0.0055 * x^2.5 / Cv^2) - CL0, [0.01, 20]);   % Eq. 3.7  mean wetted length-to-beam ratio
    CLd  = 0.012 * sqrt(lam) * tau^1.1;                                  % Eq. 3.9  C_Ld, dynamic component of the lift coefficient (first term of Eq. 3.7; Savitsky Eq. 20), used in Eq. 3.8
    Vm   = V * sqrt(max(1 - (CLd - 0.0065 * beta * CLd^0.6) / (lam * cosd(tau)), 0.36));     % Eq. 3.8  mean bottom velocity (floored at sqrt(0.36)*V)
    Re   = Vm * lam * b / c.nu;                                          % Eq. 3.10  Reynolds number
    Cf   = schoenherr(Re);                                               % Eq. 3.11 Schoenherr friction line
    Df   = (Cf + dCf) * c.rho * Vm^2 * lam * b^2 / (2 * cosd(beta));     % Eq. 3.12 (C_F,tot = C_F + dCA) inserted in Eq. 3.13 (friction drag)
    D    = W * tand(tau) + Df / cosd(tau);                               % Eq. 3.14 total hydrodynamic resistance
    Cp   = 0.75 - 1 / (5.21 * Cv^2 / lam^2 + 2.39);                      % Eq. 3.15 centre of pressure
    lc   = c.LCG - Cp * lam * b;                                         % Eq. 3.16 moment arm c (named lc to avoid clashing with the struct c)
    a    = hull.VCG - b / 4 * tand(beta);                                % Eq. 3.1  friction-line arm a
    Mres = W * ((1 - sind(tau) * sind(tau + hull.eps)) * lc / cosd(tau) - hull.f * sind(tau)) ...
           + Df * (a - hull.f);                                          % Eq. 3.17 pitching moment equilibrium (zero at tau_eq)
end
function Cf = schoenherr(Re)
    % Eq. 3.11: 0.242/sqrt(Cf) = log10(Re*Cf), solved numerically for Cf
    Cf = fzero(@(x) 0.242 / sqrt(x) - log10(Re * x), [1e-4, 0.1]);
end