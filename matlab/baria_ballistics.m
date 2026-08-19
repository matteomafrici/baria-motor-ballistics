%% ========================================================================
%  BARIA motor - complete ballistic analysis
%  BC Method (Bayern-Chemie) | Vieille's Law | c* | 0D Model | Monte Carlo
%
%  Run from the repository root:
%    matlab -batch "run('matlab/baria_ballistics.m')"
%
%  Requires: a dataset .mat file in data/ (see data/README.md), uncertainty.m
%
%  References:
%    [1] BC-V01 - Bayern-Chemie method (course document)
%    [2] Fry et al., AIAA 2001-3948 - NATO burning rate methods
%    [3] Terzic et al., Problems of Mechatronics 4(6) 2011 - 0D model
%    [4] Kallmeyer & Sayer, AIAA-82-1094 - exp vs predicted differences
%
%  See report/baria-report.pdf for the full methodology and results.
%% ========================================================================
clear; close all; clc;

%% ---- REPO ROOT (robust to run() cwd handling) ---------------------------
here = pwd;
if endsWith(here, filesep + "matlab")
    root = fileparts(here);
else
    root = here;
end
addpath(fullfile(root, 'matlab'));

%% ---- PLOT DEFAULTS: light theme -----------------------------------------
set(groot, ...
    'DefaultAxesFontSize',        11, ...
    'DefaultAxesFontName',        'Helvetica', ...
    'DefaultLineLineWidth',       1.8, ...
    'DefaultAxesBox',             'on', ...
    'DefaultAxesGridAlpha',       0.3, ...
    'DefaultFigureColor',         [0.97 0.97 0.97], ...
    'DefaultAxesColor',           [1.00 1.00 1.00], ...
    'DefaultAxesXColor',          [0.15 0.15 0.15], ...
    'DefaultAxesYColor',          [0.15 0.15 0.15], ...
    'DefaultAxesZColor',          [0.15 0.15 0.15], ...
    'DefaultAxesGridColor',       [0.15 0.15 0.15], ...
    'DefaultAxesTitleFontWeight', 'bold', ...
    'DefaultTextColor',           [0.15 0.15 0.15]);

%% ---- CONFIG: motor / propellant constants -------------------------------
rho_ap   = 1950;    % [kg/m^3]
rho_al   = 2700;    % [kg/m^3]
rho_htpb = 920;     % [kg/m^3]
rho_p    = 1 / (0.68/rho_ap + 0.18/rho_al + 0.14/rho_htpb);

d_out  = 0.160;     % grain outer diameter [m]
d_in_0 = 0.100;     % initial bore diameter [m]
l_0    = 0.290;     % grain length [m]
web_m  = (d_out - d_in_0) / 2;   % web thickness = 0.030 m

d_t = [28.80, 25.26, 21.81] * 1e-3;   % low, mid, high throat diameters
a_t = (pi/4) * d_t.^2;
dt  = 1e-3;                            % 1 ms sampling

thrsh_pct   = 0.05;   % 5% threshold (BC)
smooth_win  = 30;     % stability window [samples = ms]
dpdt_var_th = 0.12;   % normalised local std(|dp/dt|) threshold

n_mc = 5000;          % Monte Carlo draws
mc_seed = 42;         % fixed seed for reproducibility

ab_fun  = @(y) pi*(d_in_0 + 2*y)*(l_0 - 2*y) + 2*(pi/4)*(d_out^2 - (d_in_0 + 2*y)^2);
a_to_si = @(a_val, n_val) (a_val*1e-3) / (1e5)^n_val;

%% ---- LOAD DATASET (auto-detect batches in data/*.mat) -------------------
mat_files = dir(fullfile(root, 'data', '*.mat'));
if isempty(mat_files)
    error('No .mat dataset found in data/. See data/README.md.');
end

d = load(fullfile(mat_files(1).folder, mat_files(1).name));
names = sort(fieldnames(d));

pmat = {};
batch_ids = {};
num_batches = 0;
for i = 1:numel(names)
    v = d.(names{i});
    if isnumeric(v) && size(v, 2) == 3
        num_batches = num_batches + 1;
        pmat{num_batches}    = v;
        batch_ids{num_batches} = names{i};
    end
end

if num_batches == 0
    error('No N x 3 batch matrices found in the dataset file.');
end
n_traces = num_batches * 3;

%% ---- PRE-PROCESS --------------------------------------------------------
% Dataset-specific correction: restore the mid/high columns of pbar2444 and
% fix a single outlying sample (details in the report).
idx2444 = find(strcmp(batch_ids, 'pbar2444'), 1);
if ~isempty(idx2444)
    pmat{idx2444}(:, [2 3]) = pmat{idx2444}(:, [3 2]);
    pmat{idx2444}(10, 2)    = 10;
end

%% ---- COLOUR PALETTE -----------------------------------------------------
color_low = [0.82 0.10 0.10];
color_mid = [0.10 0.52 0.10];
color_hig = [0.10 0.10 0.82];
cols_regime = {color_low, color_mid, color_hig};
regimes     = {'Low','Mid','High'};
regime_names = {'low','mid','high'};

%% ---- FIGURE 1: all traces (top) + per-batch subplots (3x3) --------------
hF1 = figure('Name','BARIA - Pressure Traces', 'NumberTitle','off', ...
    'Position',[40 40 1400 980]);

ax_all = subplot(4,3,[1 2 3]);
hold(ax_all,'on'); grid(ax_all,'on');
plot(ax_all, NaN,NaN,'-','Color',color_low,'LineWidth',2,'DisplayName','Low P');
plot(ax_all, NaN,NaN,'-','Color',color_mid,'LineWidth',2,'DisplayName','Mid P');
plot(ax_all, NaN,NaN,'-','Color',color_hig,'LineWidth',2,'DisplayName','High P');
for i = 1:num_batches
    p = pmat{i};
    tms = 0:size(p,1)-1;
    for j = 1:3
        plot(ax_all, tms, p(:,j),'Color',cols_regime{j},'LineWidth',0.7,'HandleVisibility','off');
    end
end
xlabel(ax_all,'Time [ms]'); ylabel(ax_all,'Pressure [bar]');
title(ax_all,'All Experimental Pressure Traces');
legend(ax_all,'Location','northeast'); xlim(ax_all,[0 5000]);

%% ---- BC METHOD ----------------------------------------------------------
p_eff_all  = zeros(n_traces, 1);
t_burn_all = zeros(n_traces, 1);
rb_all     = zeros(n_traces, 1);
c_star_exp = zeros(n_traces, 1);
bc = struct('idx_a',{}, 'idx_g',{}, 'idx_b',{}, 'idx_e',{}, ...
            'p_ref',{}, 't_a',{}, 't_g',{}, 't_b',{}, 't_e',{});

ax_bc = gobjects(num_batches,1);
for i = 1:num_batches
    ax_bc(i) = subplot(4,3, i+3);
    hold(ax_bc(i),'on'); grid(ax_bc(i),'on');
end

count = 1;
for i = 1:num_batches
    p      = pmat{i};
    n_samp = size(p, 1);
    t      = (0:n_samp-1)' * dt;
    ax     = ax_bc(i);

    for j = 1:3
        pt = p(:, j);

        % Peak detection with 0.5 s guard window against ignition spikes.
        [pmax_j, idx_pk] = max(pt);
        if idx_pk <= 500
            [pmax_j, idx_rel] = max(pt(500:end));
            idx_pk = idx_rel + 499;
        end

        % Derivative-stability guard.
        pt_s     = movmean(pt, smooth_win);
        dpdt_s   = gradient(pt_s, dt);
        adpdt    = abs(dpdt_s);
        dpdt_var = movstd(adpdt, smooth_win) ./ (max(adpdt) + 1e-9);
        stable   = dpdt_var < dpdt_var_th;

        thrsh = thrsh_pct * pmax_j;

        % Action interval A/G.
        cands_rise = find(pt(1:idx_pk) >= thrsh);
        if isempty(cands_rise)
            idx_a = 1;
        else
            c0 = cands_rise(1);
            sf = find(stable(c0:idx_pk)) + c0 - 1;
            if ~isempty(sf), idx_a = sf(1); else, idx_a = c0; end
        end

        cands_fall = find(pt(idx_pk:end) >= thrsh) + idx_pk - 1;
        if isempty(cands_fall)
            idx_g = n_samp;
        else
            c1 = cands_fall(end);
            sb = find(stable(idx_pk:c1)) + idx_pk - 1;
            if ~isempty(sb), idx_g = sb(end); else, idx_g = c1; end
        end

        t_a = t(idx_a); t_g = t(idx_g);

        % Reference pressure (BC).
        i1    = trapz(t(idx_a:idx_g), pt(idx_a:idx_g)) / 2;
        p_ref = i1 / (t_g - t_a);

        % Burning interval B/E: maximal region with p >= p_ref in [A,G].
        mask_pref = false(n_samp, 1);
        mask_pref(idx_a:idx_g) = pt(idx_a:idx_g) >= p_ref;
        idx_be = find(mask_pref);
        if isempty(idx_be)
            idx_b = idx_a; idx_e = idx_g;
        else
            idx_b = idx_be(1); idx_e = idx_be(end);
        end
        t_b = t(idx_b); t_e = t(idx_e);
        t_burn_j = t_e - t_b;

        p_eff_j = trapz(t(idx_b:idx_e), pt(idx_b:idx_e)) / t_burn_j;   % [bar]
        rb_j    = web_m * 1e3 / t_burn_j;                              % [mm/s]
        c_star_exp(count) = (p_eff_j*1e5 * t_burn_j * a_t(j)) / ...
            (rho_p * (pi/4) * (d_out^2 - d_in_0^2) * l_0);             % [m/s]

        p_eff_all(count)  = p_eff_j;
        t_burn_all(count) = t_burn_j;
        rb_all(count)     = rb_j;
        bc(count).idx_a = idx_a; bc(count).idx_g = idx_g;
        bc(count).idx_b = idx_b; bc(count).idx_e = idx_e;
        bc(count).p_ref = p_ref;
        bc(count).t_a = t_a; bc(count).t_g = t_g;
        bc(count).t_b = t_b; bc(count).t_e = t_e;

        plot(ax, t*1e3, pt,'Color',cols_regime{j},'LineWidth',1.4, ...
            'DisplayName',sprintf('%s  P_{eff}=%.2f bar', regime_names{j}, p_eff_j));
        xline(ax, t_a*1e3,'--','Color',cols_regime{j},'Alpha',0.5,'HandleVisibility','off');
        xline(ax, t_g*1e3,'--','Color',cols_regime{j},'Alpha',0.5,'HandleVisibility','off');
        xline(ax, t_b*1e3,'-.','Color',cols_regime{j},'Alpha',0.5,'HandleVisibility','off');
        xline(ax, t_e*1e3,'-.','Color',cols_regime{j},'Alpha',0.5,'HandleVisibility','off');
        yline(ax, p_ref,   ':','Color',cols_regime{j},'Alpha',0.5,'HandleVisibility','off');

        count = count + 1;
    end
    legend(ax,'Location','northeast','FontSize',7);
    xlabel(ax,'Time [ms]'); ylabel(ax,'P [bar]');
    title(ax, batch_ids{i});
end

exportgraphics(hF1, fullfile(root, 'figures', 'pressure_traces.png'), ...
    'Resolution', 200, 'BackgroundColor', 'white');

%% ---- PER-REGIME STATISTICS ---------------------------------------------
fprintf('\n========== RESULTS PER PRESSURE REGIME ==========\n');
for j = 1:3
    idx_reg = j:3:n_traces;
    P = p_eff_all(idx_reg);  R = rb_all(idx_reg);  T = t_burn_all(idx_reg);
    fprintf('--- %s Pressure (D_t = %.2f mm) ---\n', regimes{j}, d_t(j)*1e3);
    fprintf('  P_eff  : mean = %6.3f bar   CV = %.2f%%\n', mean(P), 100*std(P)/mean(P));
    fprintf('  r_b    : mean = %6.4f mm/s  CV = %.2f%%\n', mean(R), 100*std(R)/mean(R));
    fprintf('  t_burn : mean = %6.4f s     sigma = %5.4f s\n', mean(T), std(T));
end

%% ---- VIEILLE'S LAW FIT & c* --------------------------------------------
[a, inc_a, n, inc_n, r2] = uncertainty(p_eff_all, rb_all);
c_star_mean = mean(c_star_exp);
c_star_std  = std(c_star_exp);

fprintf('\n========== VIEILLE''S LAW ==========\n');
fprintf('  a = %.5f +/- %.5f [mm/s/bar^n]\n', a, inc_a);
fprintf('  n = %.5f +/- %.5f\n', n, inc_n);
fprintf('  R2 = %.6f\n', r2);
fprintf('\n========== CHARACTERISTIC VELOCITY ==========\n');
fprintf('  c* = %.2f +/- %.2f m/s  (CV = %.2f%%)\n', ...
    c_star_mean, c_star_std, 100*c_star_std/c_star_mean);

%% ---- FIGURE 2: Vieille's law (log-log, MC +/-1 sigma band) -------------
hF2 = figure('Name','Vieille Law - Report', 'NumberTitle','off', ...
    'Position',[200 200 620 450]);

ax_vr = subplot(1,1,1);
hold(ax_vr,'on'); grid(ax_vr,'on');

p_fit_r  = logspace(log10(min(p_eff_all)-2), log10(max(p_eff_all)+2), 400);
rb_fit_r = a         .* p_fit_r.^n;
rb_up_r  = (a+inc_a) .* p_fit_r.^(n+inc_n);
rb_dn_r  = (a-inc_a) .* p_fit_r.^(n-inc_n);

fill(ax_vr, [p_fit_r, fliplr(p_fit_r)], [rb_up_r, fliplr(rb_dn_r)], ...
    [0.75 0.75 0.75], 'FaceAlpha', 0.45, 'EdgeColor','none', ...
    'DisplayName', 'MC $\pm1\sigma$ band');

scatter(ax_vr, p_eff_all(1:3:end), rb_all(1:3:end), 55, color_low, 'filled', ...
    'DisplayName', sprintf('Low P ($\\bar{P}_{eff}=%.1f$ bar)', mean(p_eff_all(1:3:end))));
scatter(ax_vr, p_eff_all(2:3:end), rb_all(2:3:end), 55, color_mid, 'filled', ...
    'DisplayName', sprintf('Mid P ($\\bar{P}_{eff}=%.1f$ bar)', mean(p_eff_all(2:3:end))));
scatter(ax_vr, p_eff_all(3:3:end), rb_all(3:3:end), 55, color_hig, 'filled', ...
    'DisplayName', sprintf('High P ($\\bar{P}_{eff}=%.1f$ bar)', mean(p_eff_all(3:3:end))));

plot(ax_vr, p_fit_r, rb_fit_r, 'k-', 'LineWidth', 2.2, ...
    'DisplayName', sprintf('$r_b = %.4f\\,P^{%.4f}$  ($R^2=%.4f$)', a, n, r2));

set(ax_vr, 'XScale','log', 'YScale','log');
xlabel(ax_vr, '$P_{\mathrm{eff}}$ [bar]', 'Interpreter','latex');
ylabel(ax_vr, '$r_b$ [mm/s]',             'Interpreter','latex');
title(ax_vr,  "Vieille's Law - AP/HTPB/Al propellant");
legend(ax_vr, 'Location','northwest', 'Interpreter','latex');
set(ax_vr, 'XTick',[30 40 50 60 70 80], 'XTickLabel',{'30','40','50','60','70','80'});
set(ax_vr, 'YTick',[6 7 8 9 10],        'YTickLabel',{'6','7','8','9','10'});
xlim(ax_vr, [28 82]);
ylim(ax_vr, [6.0 9.8]);

exportgraphics(hF2, fullfile(root, 'figures', 'vieille_fit.png'), ...
    'Resolution', 200, 'BackgroundColor', 'white');

%% ---- FIGURE 3: 0D model, all throats + bounds on a ----------------------
hF3 = figure('Name','0D Model - All throats + bounds', 'NumberTitle','off', ...
    'Position',[100 100 900 500]);
hold on; grid on;

a_cases = [a - inc_a, a, a + inc_a];
line_styles = {'--','-','--'};
lw_cases    = [1.2, 2.2, 1.2];

for j = 1:3
    for k = 1:3
        a_si = a_to_si(a_cases(k), n);
        y = 0; t_now = 0; t_sim = []; p_sim = [];
        while y < web_m
            lc = l_0 - 2*y;
            if lc <= 0, break; end
            ab   = ab_fun(y);
            p_pa = (a_si * rho_p * c_star_mean * (ab/a_t(j)))^(1/(1-n));
            t_sim(end+1) = t_now;
            p_sim(end+1) = p_pa / 1e5;
            y     = y     + a_si * p_pa^n * dt;
            t_now = t_now + dt;
        end
        hv = 'off';
        if k == 2, hv = 'on'; end
        plot(t_sim, p_sim, line_styles{k}, 'Color', cols_regime{j}, ...
            'LineWidth', lw_cases(k), 'HandleVisibility', hv, ...
            'DisplayName', sprintf('D_t=%.2f mm', d_t(j)*1e3));
    end
    plot(NaN, NaN, '-', 'Color', cols_regime{j}, 'LineWidth', 2.2, ...
        'DisplayName', sprintf('D_t=%.2f mm (nom. + bounds)', d_t(j)*1e3));
end
xlabel('Time [s]'); ylabel('Chamber Pressure [bar]');
title('0D Quasi-Steady Model - all nozzle configs (bounds on a)');
legend('Location','northeast');

exportgraphics(hF3, fullfile(root, 'figures', 'ballistic_model.png'), ...
    'Resolution', 200, 'BackgroundColor', 'white');

%% ---- EXP vs MODEL + RESIDUALS (first batch, mid pressure) ---------------
batch_idx = 1;
col_idx   = 2;
p_exp_raw = pmat{batch_idx}(:, col_idx);
t_exp     = (0:length(p_exp_raw)-1)' * dt;

a_si   = a_to_si(a, n);
at_med = a_t(col_idx);
y = 0; t_now = 0; t_sim = []; p_sim = [];
while y < web_m
    lc = l_0 - 2*y;
    if lc <= 0, break; end
    ab   = ab_fun(y);
    p_pa = (a_si * rho_p * c_star_mean * (ab/at_med))^(1/(1-n));
    t_sim(end+1) = t_now;
    p_sim(end+1) = p_pa / 1e5;
    y     = y     + a_si * p_pa^n * dt;
    t_now = t_now + dt;
end

p_model_interp = interp1(t_sim, p_sim, t_exp, 'linear', 0);
residual       = p_exp_raw - p_model_interp;

c_med = (batch_idx-1)*3 + col_idx;
t_win = t_exp >= bc(c_med).t_b & t_exp <= bc(c_med).t_e;
rmse  = sqrt(mean(residual(t_win).^2));
fprintf('\nExp vs Model (%s, mid P): RMSE_{B->E} = %.3f bar\n', batch_ids{batch_idx}, rmse);

hF4 = figure('Name','Exp vs Model + Residuals', 'NumberTitle','off', ...
    'Position',[100 100 760 620]);

subplot(3,1,[1 2]);
hold on; grid on;
plot(t_exp*1e3, p_exp_raw, 'Color', color_mid, 'LineWidth', 1.2, 'DisplayName','Experimental');
plot(t_sim*1e3, p_sim, 'r--', 'LineWidth', 2.0, 'DisplayName','0D Model (nominal)');
xline(bc(c_med).t_b*1e3, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
xline(bc(c_med).t_e*1e3, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
ylabel('Pressure [bar]');
legend('Location','south');
title(sprintf('%s - Mid pressure: Experimental vs 0D Model', batch_ids{batch_idx}));

subplot(3,1,3);
hold on; grid on;
plot(t_exp*1e3, residual, 'k-', 'LineWidth', 1);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Time [ms]'); ylabel('\DeltaP [bar]');
title(sprintf('Residuals (exp - model) | RMSE_{B->E} = %.3f bar', rmse));

exportgraphics(hF4, fullfile(root, 'figures', 'exp_vs_model.png'), ...
    'Resolution', 200, 'BackgroundColor', 'white');

%% ---- MONTE CARLO --------------------------------------------------------
fprintf('\n========== MONTE CARLO (N=%d, seed %d) ==========\n', n_mc, mc_seed);
rng(mc_seed);

% a and n: aleatoric; c*: epistemic.
a_mc     = a           + inc_a      * randn(n_mc, 1);
n_mc_vec = n           + inc_n      * randn(n_mc, 1);
cstar_mc = c_star_mean + c_star_std * randn(n_mc, 1);

% Independent shuffle to break sampling-order correlation.
a_mc     = a_mc(randperm(n_mc));
n_mc_vec = n_mc_vec(randperm(n_mc));
cstar_mc = cstar_mc(randperm(n_mc));

mc = struct();

hF5 = figure('Name','Monte Carlo', 'NumberTitle','off', ...
    'Position',[100 20 1150 700]);
ax_cv = subplot(2,3,[1 2 3]);
hold(ax_cv,'on'); grid(ax_cv,'on');

for jj = 1:3
    t_burn_mc = zeros(n_mc, 1);
    p_eff_mc  = zeros(n_mc, 1);
    t_conv    = zeros(n_mc, 1);

    for ii = 1:n_mc
        a_si_ii = a_to_si(a_mc(ii), n_mc_vec(ii));

        y = 0;   t_now = 0;   sum_p_dt = 0;
        while y < web_m
            lc = l_0 - 2*y;
            if lc <= 0, break; end
            ab    = ab_fun(y);
            p_pa  = (a_si_ii * rho_p * cstar_mc(ii) * (ab/a_t(jj)))^(1/(1 - n_mc_vec(ii)));
            rb_si = a_si_ii * p_pa^n_mc_vec(ii);
            sum_p_dt = sum_p_dt + (p_pa/1e5) * dt;
            y     = y     + rb_si * dt;
            t_now = t_now + dt;
        end

        t_burn_mc(ii) = t_now;
        p_eff_mc(ii)  = sum_p_dt / t_now;
        t_conv(ii)    = mean(t_burn_mc(1:ii));
    end

    mc(jj).t_burn  = t_burn_mc;
    mc(jj).p_eff   = p_eff_mc;
    mc(jj).mean_tb = mean(t_burn_mc);
    mc(jj).std_tb  = std(t_burn_mc);
    mc(jj).cv_tb   = mc(jj).std_tb / mc(jj).mean_tb * 100;
    mc(jj).mean_pe = mean(p_eff_mc);
    mc(jj).cv_pe   = std(p_eff_mc) / mc(jj).mean_pe * 100;

    fprintf('--- %s P  D_t=%.2f mm ---\n', regimes{jj}, d_t(jj)*1e3);
    fprintf('  t_burn : mean = %.4f s  sigma = %.4f s  CV = %.2f%%\n', ...
        mc(jj).mean_tb, mc(jj).std_tb, mc(jj).cv_tb);
    fprintf('  P_eff  : mean = %.3f bar  CV = %.2f%%\n', ...
        mc(jj).mean_pe, mc(jj).cv_pe);

    plot(ax_cv, 1:n_mc, t_conv, 'Color', cols_regime{jj}, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('D_t=%.2f mm', d_t(jj)*1e3));
end

xlabel(ax_cv,'Iterations'); ylabel(ax_cv,'mean t_b [s]');
title(ax_cv,'MC convergence - running mean of t_{burn}');
legend(ax_cv,'Location','east');

for jj = 1:3
    ax_h = subplot(2,3, 3+jj);
    hold(ax_h,'on'); grid(ax_h,'on');
    histogram(ax_h, mc(jj).t_burn, 60, 'FaceColor', cols_regime{jj}, ...
        'EdgeColor','w', 'FaceAlpha', 0.85);
    xline(ax_h, mc(jj).mean_tb,               'k-',  'LineWidth', 2.0);
    xline(ax_h, mc(jj).mean_tb + mc(jj).std_tb, 'k--', 'LineWidth', 1.2);
    xline(ax_h, mc(jj).mean_tb - mc(jj).std_tb, 'k--', 'LineWidth', 1.2);
    xlabel(ax_h,'t_b [s]'); ylabel(ax_h,'Count');
    title(ax_h, sprintf('%s P  D_t=%.2f mm\n\\mu=%.4f s   \\sigma=%.4f s   CV=%.2f%%', ...
        regimes{jj}, d_t(jj)*1e3, mc(jj).mean_tb, mc(jj).std_tb, mc(jj).cv_tb));
end

exportgraphics(hF5, fullfile(root, 'figures', 'monte_carlo.png'), ...
    'Resolution', 200, 'BackgroundColor', 'white');

%% ---- FINAL SUMMARY ------------------------------------------------------
fprintf('\n========== FINAL SUMMARY ==========\n');
fprintf('rho_p  = %.1f kg/m3\n', rho_p);
fprintf('m_tot  = %.5f kg\n',    rho_p * (pi/4) * (d_out^2 - d_in_0^2) * l_0);
fprintf('web    = %.1f mm\n',    web_m * 1e3);
fprintf('a      = %.5f +/- %.5f [mm/s/bar^n]\n', a, inc_a);
fprintf('n      = %.5f +/- %.5f\n', n, inc_n);
fprintf('R2     = %.6f\n', r2);
fprintf('c*     = %.2f +/- %.2f m/s  (CV=%.2f%%)\n', ...
    c_star_mean, c_star_std, 100*c_star_std/c_star_mean);