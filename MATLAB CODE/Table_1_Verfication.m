% Table_1_Verfication.m
% Checks my Savitsky (1964) code against the worked example in
% Savitsky (1964), Table 1, p. 89 (Chapter 3, Table 3.1).
% Units: Imperial, as in the paper (lb, ft, s, slug).
clear; clc;
% Input data from Savitsky (1964), Table 1
b    = 14;        % chine beam (ft)
W    = 60000;     % displacement weight (lb)
beta = 10;        % deadrise angle (deg)
LCG  = 29.0;      % LCG from transom (ft)
VCG  = 2.0;       % VCG above keel (ft)
f    = 0.50;      % thrust line to CG (ft)
epsilon = 4;      % thrust inclination to keel (deg)
V    = 67.5;      % speed, 40 knots (ft/s)
g    = 32.2;      % gravitational acceleration (ft/s^2)
rho  = 1.94;      % water density (slug/ft^3)
nu   = 1.0e-5;    % kinematic viscosity (ft^2/s)
dCf  = 0.0004;    % ATTC roughness allowance, Delta C_f (Eq. 3.12)
%  Values printed in Savitsky's Table 1 at trim = 2, 3 and 4 deg 
tauTrial = [2 3 4];
paper.lambda = [3.85     2.60     1.86   ];   % Row 4
paper.Vm     = [67.0     66.6     66.2   ];   % Row 5
paper.Re     = [3.61e8   2.42e8   1.73e8 ];   % Row 6
paper.Cf     = [0.00174  0.00184  0.00192];   % Row 7
paper.Df     = [7340     5160     3760   ];   % Row 10
paper.D      = [9434     8304     7948   ];   % Row 16
paper.Mres   = [-149960  336600   648350 ];   % Row 30 (moment residual)
%  Run my code at the same three trim angles 
for i = 1:3
    r = planing(tauTrial(i), b, W, beta, LCG, VCG, f, epsilon, V, g, rho, nu, dCf);
    mine.lambda(i) = r.lambda;
    mine.Vm(i)     = r.Vm;
    mine.Re(i)     = r.Re;
    mine.Cf(i)     = r.Cf;
    mine.Df(i)     = r.Df;
    mine.D(i)      = r.D;
    mine.Mres(i)   = r.Mres;
end
% Print the row-by-row comparison 
rows   = {'lambda', 'Vm', 'Re', 'Cf', 'Df', 'D', 'Mres'};
labels = {'lambda (Eq 3.7)', 'Vm (Eq 3.8)', 'Re (Eq 3.10)', 'Cf (Eq 3.11)', ...
          'Df (Eq 3.13)', 'D (Eq 3.14)', 'Mres (Eq 3.17)'};
fprintf('%-16s', 'Quantity');
for i = 1:3
    fprintf('|  tau = %d deg: mine, paper, diff   ', tauTrial(i));
end
fprintf('\n');
for k = 1:numel(rows)
    fprintf('%-16s', labels{k});
    for i = 1:3
        m = mine.(rows{k})(i);
        p = paper.(rows{k})(i);
        fprintf('|  %10.4g %10.4g %+6.1f%%  ', m, p, percentDiff(m, p));
    end
    fprintf('\n');
end
%    Equilibrium trim and drag (Eq. 3.18) 
% Linear interpolation between 2 and 3 deg, where the moment residual
tauEq = interpTrim(tauTrial(1), tauTrial(2), mine.Mres(1), mine.Mres(2));
DEq   = mine.D(1) + (tauEq - tauTrial(1)) * (mine.D(2) - mine.D(1));

% Savitsky's equilibrium trim and drag, as printed in Table 1
tauPaper = 2.3;                                              % deg
DPaper   = 9095;                                             % lb

fprintf('\nEquilibrium (linear interpolation between 2 and 3 deg, Eq. 3.18):\n');
fprintf('  trim: mine %.3f deg, paper %.1f deg, diff %+.1f%%\n', ...
        tauEq, tauPaper, percentDiff(tauEq, tauPaper));
fprintf('  drag: mine %.0f lb,   paper %.0f lb,   diff %+.1f%%\n', ...
        DEq, DPaper, percentDiff(DEq, DPaper));
% Exact equilibrium (no interpolation)
residual = @(t) getResidual(t, b, W, beta, LCG, VCG, f, epsilon, V, g, rho, nu, dCf);
tauExact = fzero(residual, [2 3]);
rExact   = planing(tauExact, b, W, beta, LCG, VCG, f, epsilon, V, g, rho, nu, dCf);
fprintf('\nExact equilibrium (root of Eq. 3.17): trim %.3f deg, drag %.0f lb\n', ...
        tauExact, rExact.D);
% Side check: Savitsky's Row 10 from his own numbers 
% Df at 2 deg from Eq. 3.13, using his Rows 4, 5 and 9 (Cf + dCf = 0.00214)
Df10 = 0.00214 * rho * 67.0^2 * 3.85 * b^2 / (2 * cosd(beta));
fprintf('\nRow 10 at 2 deg recomputed from the paper''s own values: %.0f lb (printed 7340 lb)\n', Df10);
% Functions
function r = planing(tau, b, W, beta, LCG, VCG, f, epsilon, V, g, rho, nu, dCf)
% Savitsky (1964) planing equations at one trial trim angle tau (deg).
% Returns a struct with the main results and the moment residual Mres
% (Eq. 3.17), which is zero at the equilibrium trim.
    % Eq. 3.4: speed coefficient
    Cv = V / sqrt(g * b);
    % Eq. 3.5: lift coefficient of the deadrise surface
    CLb = W / (0.5 * rho * V^2 * b^2);
    % Eq. 3.6: zero-deadrise lift coefficient (solve CLb = CL0 - 0.0065*beta*CL0^0.6)
    CL0 = fzero(@(x) x - 0.0065 * beta * x^0.6 - CLb, [CLb, 2 * CLb + 0.1]);
    % Eq. 3.7: mean wetted length-to-beam ratio lambda
    lambda = fzero(@(x) tau^1.1 * (0.012 * sqrt(x) + 0.0055 * x^2.5 / Cv^2) - CL0, ...
                   [0.01, 20]);
    % Eq. 3.9: dynamic part of the lift coefficient
    CLd = 0.012 * sqrt(lambda) * tau^1.1;
    % Eq. 3.8: mean bottom velocity (term under the root floored at 0.36)
    term = 1 - (CLd - 0.0065 * beta * CLd^0.6) / (lambda * cosd(tau));
    Vm   = V * sqrt(max(term, 0.36));
    % Eq. 3.10: Reynolds number based on the mean wetted length
    Re = Vm * lambda * b / nu;
    % Eq. 3.11: Schoenherr friction coefficient
    Cf = schoenherr(Re);
    % Eqs. 3.12 and 3.13: friction drag, with the roughness allowance added
    CfTot = Cf + dCf;
    Df    = CfTot * rho * Vm^2 * lambda * b^2 / (2 * cosd(beta));
    % Eq. 3.14: total resistance
    D = W * tand(tau) + Df / cosd(tau);
    % Eq. 3.15: centre of pressure, as a fraction of the mean wetted length
    Cp = 0.75 - 1 / (5.21 * Cv^2 / lambda^2 + 2.39);
    % Eq. 3.16: moment arm c between the centre of pressure and the CG
    c = LCG - Cp * lambda * b;
    % Eq. 3.1: distance a between the friction force and the CG
    a = VCG - (b / 4) * tand(beta);
    % Eq. 3.17: pitching-moment residual about the CG
    Mres = W * ((1 - sind(tau) * sind(tau + epsilon)) * c / cosd(tau) - f * sind(tau)) ...
           + Df * (a - f);
    % Collect the results
    r.lambda = lambda;
    r.Vm     = Vm;
    r.Re     = Re;
    r.Cf     = Cf;
    r.Df     = Df;
    r.D      = D;
    r.Mres   = Mres;
end
function Cf = schoenherr(Re)
% Eq. 3.11: 0.242/sqrt(Cf) = log10(Re*Cf), solved for Cf
    Cf = fzero(@(x) 0.242 / sqrt(x) - log10(Re * x), [1e-4, 0.1]);
end
function M = getResidual(tau, b, W, beta, LCG, VCG, f, epsilon, V, g, rho, nu, dCf)
% Moment residual only, so that fzero can find the equilibrium trim
    r = planing(tau, b, W, beta, LCG, VCG, f, epsilon, V, g, rho, nu, dCf);
    M = r.Mres;
end
function tauEq = interpTrim(tau1, tau2, M1, M2)
% Eq. 3.18: linear interpolation of the trim where the residual is zero
    tauEq = tau1 - M1 * (tau2 - tau1) / (M2 - M1);
end
function d = percentDiff(mine, paper)
% Percentage difference of my value relative to the paper's value
    d = 100 * (mine - paper) / abs(paper);
end