%% ------------------------------------------------------------------
%  variogram_fit_seasonal.m
%  Part of: radar-gauge-precipitation (Payam Heidarian, University of
%  Brescia / MeteoSwiss collaboration). Inline comments are partly in
%  Persian. Data are not included; set the <DATA_ROOT> paths below to
%  your own folders before running.
%% ------------------------------------------------------------------
%% ============================================================
%  DIAGNOSTIC PART B -- PER SEASON (full pipeline, repeated x4)
%  ------------------------------------------------------------
%  This is the EXACT SAME pipeline as the uploaded annual
%  "Diagnostic Part B" script (same formulas, same 7 figures,
%  same 2 tables, same custom tooltips) -- the ONLY difference is
%  that the pair-computation step (PART 3) is restricted to one
%  season's hours at a time, and outputs go to a season-specific
%  subfolder. Nothing about the estimation logic, the WLS fitting,
%  or the figure/table contents has been changed or simplified.
%% ============================================================
clc; clear;

USE_BINS = false;   % kept identical to the uploaded script's setting
N_BINS   = 705;

path_gs = "<DATA_ROOT>\" + ...
    "Profession\1-PRECIPITATION\Code\CODE\Moghayese gs vs nc franciacorta\" + ...
    "چیزی که استاد خواسته\code\1-Combi-GS-Valid H\1-in\GS-cel-xy-val H-pre";

csv6_names = {
    'Consorzio_Clusane_Belvedere_6620_corretto'
    'Consorzio_Sant_Eusebio_40485_corretto'
    'Consorzio_Brescia_Sant_Anna_6619_corretto'
    'Consorzio_Dosso_Baiana_40484_corretto'
    'Consorzio_Crociate_40479_corretto'
    'Consorzio_Versella_40480_corretto'
};

season_names = {'Spring','Summer','Autumn','Winter'};
season_map   = [4,4,1,1,1,2,2,2,3,3,3,4];   % month -> season index (1..4)

%% =====================================================================
%% PART 1 -- read data (identical to the uploaded script)
%% =====================================================================
all_files   = dir(fullfile(path_gs, 'Consorzio*.csv'));
all_names   = erase({all_files.name}, '.csv');
train_names = setdiff(all_names, csv6_names, 'stable');
n_files     = numel(train_names);

raw = cell(n_files, 1);
cnt = 0;
for i = 1:n_files
    fpath = fullfile(path_gs, [train_names{i} '.csv']);
    if ~isfile(fpath), continue; end
    T = readtable(fpath);
    T = T(:, 1:7);
    T.Properties.VariableNames = {'Time','Sx','Sy','Rx','Ry','StRain','RaRain'};
    T.Time = datetime(T.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    valid  = T.StRain > 0 & T.RaRain > 0;
    if sum(valid) < 2, continue; end
    cnt = cnt + 1;
    raw{cnt}.Sx   = T.Sx(1);
    raw{cnt}.Sy   = T.Sy(1);
    raw{cnt}.Time = T.Time(valid);
    raw{cnt}.lnr  = log(T.StRain(valid) ./ T.RaRain(valid));
end
raw = raw(1:cnt);
n   = cnt;
fprintf('Stations loaded: %d\n', n);

%% =====================================================================
%% PART 2 -- dense grid (identical to the uploaded script)
%% =====================================================================
grid_start = min(cellfun(@(r) min(r.Time), raw));
grid_end   = max(cellfun(@(r) max(r.Time), raw));
grid_times = (grid_start:hours(1):grid_end)';
nGrid      = numel(grid_times);

Sx = cellfun(@(r) r.Sx, raw)';
Sy = cellfun(@(r) r.Sy, raw)';

lnr_dense = nan(nGrid, n);
for s = 1:n
    idx = round(hours(raw{s}.Time - grid_start)) + 1;
    ok  = idx >= 1 & idx <= nGrid;
    lnr_dense(idx(ok), s) = raw{s}.lnr(ok);
end

Dfull = squareform(pdist([Sx(:), Sy(:)]));
fprintf('Dense grid: %d hours\n', nGrid);

%% NEW (only addition vs. the uploaded script): a season mask per hour
season_of_grid = season_map(month(grid_times))';

%% =====================================================================
%% PART 7 (model definitions) -- identical, computed once, shared by
%% all seasons since the model FORMS themselves don't depend on season
%% =====================================================================
exp_model = @(p,h) p(1) + p(2).*(1 - exp(-h./p(3)));
sph_model = @(p,h) p(1) + p(2).* ...
    ((1.5*(h./p(3)) - 0.5*(h./p(3)).^3).*(h <= p(3)) + ...
     1.*(h > p(3)));
gau_model = @(p,h) p(1) + p(2).*(1 - exp(-(h./p(3)).^2));
model_list  = {exp_model, sph_model, gau_model};
model_names = {'Exponential', 'Spherical', 'Gaussian'};

opts = optimoptions('lsqnonlin', ...
    'Display',            'off', ...
    'MaxIterations',      5000, ...
    'FunctionTolerance',  1e-10, ...
    'StepTolerance',      1e-10, ...
    'OptimalityTolerance',1e-10);

lb = [0,     0.05,   500   ];
ub = [0.5,   4.0,    200000];

%% =====================================================================
%% MAIN LOOP -- run the ENTIRE PART 3-13 pipeline once per season
%% =====================================================================
for season_idx = 1:4
    season_mask = (season_of_grid == season_idx);
    season_label = season_names{season_idx};

    if USE_BINS
        out_dir = sprintf('Output_Diagnostic_B_Binned_%s', season_label);
    else
        out_dir = sprintf('Output_Diagnostic_B_RawPairs_%s', season_label);
    end
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    fprintf('\n\n========== SEASON: %s ==========\n', season_label);

    %% =================================================================
    %% PART 3 -- empirical variogram for all station pairs, RESTRICTED
    %% to this season's hours (the only logic change vs. the uploaded
    %% script: "mask" now also requires season_mask to be true)
    %% =================================================================
    fprintf('Computing empirical variogram for all station pairs (season-restricted)...\n');

    pair_h = [];
    pair_g = [];
    pair_n = [];

    for ii = 1:n-1
        for jj = ii+1:n
            dh   = Dfull(ii,jj);
            mask = season_mask & ~isnan(lnr_dense(:,ii)) & ~isnan(lnr_dense(:,jj));
            np   = sum(mask);
            if np < 10, continue; end
            diffs  = lnr_dense(mask,ii) - lnr_dense(mask,jj);
            g_pair = sum(diffs.^2) / (2*np);
            pair_h(end+1,1) = dh; %#ok<AGROW>
            pair_g(end+1,1) = g_pair; %#ok<AGROW>
            pair_n(end+1,1) = np; %#ok<AGROW>
        end
    end
    n_pairs = numel(pair_h);
    fprintf('Station pairs computed: %d\n', n_pairs);

    %% =================================================================
    %% PART 4 -- raw vs binned display points (identical logic)
    %% =================================================================
    if USE_BINS
        h_max  = prctile(pair_h, 95);
        edges_h = linspace(0, h_max, N_BINS+1);
        bin_h = nan(N_BINS,1); bin_g = nan(N_BINS,1); bin_n = zeros(N_BINS,1);
        for b = 1:N_BINS
            mask     = pair_h >= edges_h(b) & pair_h < edges_h(b+1);
            bin_n(b) = sum(mask);
            if bin_n(b) < 2, continue; end
            bin_h(b) = mean(pair_h(mask));
            bin_g(b) = sum(pair_n(mask) .* pair_g(mask)) / sum(pair_n(mask));
        end
        valid_bins = ~isnan(bin_h);
        plot_h = bin_h(valid_bins); plot_g = bin_g(valid_bins); plot_n = bin_n(valid_bins);
        scatter_label = 'Empirical  \gamma(h)   (size \propto n pairs,  binned)';
        mode_str      = sprintf('Binned  (%d bins)', N_BINS);
    else
        h_max = prctile(pair_h, 95);
        plot_h = pair_h; plot_g = pair_g; plot_n = pair_n;
        scatter_label = sprintf('Empirical  \\gamma(h)   (%d station pairs,  size \\propto n_{shared hours},  n_{min}=%d  n_{max}=%d)', ...
            n_pairs, min(pair_n), max(pair_n));
        mode_str      = 'Raw pairs (no binning)';
    end

    %% =================================================================
    %% PART 5 -- initial parameter estimates (identical logic)
    %% =================================================================
    short_threshold = prctile(pair_h, 10);
    short_mask      = pair_h <= short_threshold;
    C0_init         = max(mean(pair_g(short_mask)), 0.005);

    sill_total_init = prctile(pair_g, 85);
    sill_total_init = max(sill_total_init, C0_init + 0.05);
    CH_init         = sill_total_init - C0_init;
    target_63       = C0_init + 0.63 * CH_init;

    if USE_BINS
        bg_tmp = bin_g(valid_bins);
        bh_tmp = bin_h(valid_bins);
    else
        tmp_edges = linspace(0, h_max, 51);
        bh_tmp = []; bg_tmp = [];
        for b = 1:50
            mk = pair_h >= tmp_edges(b) & pair_h < tmp_edges(b+1);
            if sum(mk) < 1, continue; end
            bh_tmp(end+1) = mean(pair_h(mk)); %#ok<AGROW>
            bg_tmp(end+1) = mean(pair_g(mk)); %#ok<AGROW>
        end
    end
    [~, range_bin_idx] = min(abs(bg_tmp - target_63));
    AH_init = bh_tmp(range_bin_idx);
    AH_init = max(AH_init, 500);

    p0 = [C0_init, CH_init, AH_init];
    p0 = min(max(p0, lb+1e-6), ub-1e-6);

    fprintf('\n--- Initial parameter estimates (%s) ---\n', season_label);
    fprintf('Nugget  C0  = %.4f\n', C0_init);
    fprintf('Sill    CH  = %.4f\n', CH_init);
    fprintf('Range   AH  = %.1f km\n', AH_init/1000);

    %% =================================================================
    %% PART 6 -- Fig_B0: initial-parameter diagnostic (identical)
    %% =================================================================
    h_plot  = linspace(0, h_max, 600);
    colors3 = {[0.20 0.47 0.75], [0.85 0.33 0.10], [0.47 0.67 0.19]};

    if USE_BINS
        sz_scale = 15 + 80*(plot_n / max(plot_n));
    else
        sz_scale = 4 + 18*(plot_n / max(plot_n));
    end
    sz_display = round(sz_scale);

    fig0 = figure('Position', [50 50 860 540]);
    hold on;
    sc0 = scatter(plot_h/1000, plot_g, sz_display, 'k', 'filled', 'DisplayName', scatter_label);
    sc0.UserData = plot_n;
    if ~USE_BINS, sc0.MarkerFaceAlpha = 0.4; sc0.MarkerEdgeAlpha = 0; end

    yline(C0_init, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.2, ...
        'DisplayName', sprintf('Nugget  C_0 = %.4f', C0_init));
    yline(C0_init + CH_init, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 2.2, ...
        'DisplayName', sprintf('Total sill  C_0+C_H = %.4f', C0_init+CH_init));
    yline(target_63, ':', 'Color', [0.20 0.47 0.75], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(AH_init/1000, '-', 'Color', [0.20 0.47 0.75], 'LineWidth', 2.2, ...
        'DisplayName', sprintf('Range  A_H \\approx %.1f km', AH_init/1000));
    xline(short_threshold/1000, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, ...
        'DisplayName', sprintf('10th pct of distances = %.1f km', short_threshold/1000));

    y_range = prctile(plot_g, 97);

    annotation('textbox', [0.52 0.18 0.40 0.28], ...
        'String', {
            'Parameter estimation (Chiles & Delfiner, 2012):', ...
            sprintf('  C_0 = mean[\\gamma(h)] for h < %.1f km', short_threshold/1000), ...
            '  C_0 = nugget effect (micro-scale variability)', ...
            ' ', ...
            '  C_0+C_H = 85th pct of \\gamma(h)  [Sill]', ...
            '  C_H = partial sill', ...
            ' ', ...
            '  A_H: h where \\gamma(A_H) = C_0 + 0.632·C_H', ...
            sprintf('  Mode: %s', mode_str), ...
            sprintf('  Season: %s', season_label)}, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'w', ...
        'EdgeColor', [0.5 0.5 0.5], 'FontSize', 8.5);

    xlabel('Lag distance  h  (km)', 'FontSize', 12);
    ylabel('\gamma(h)  --  variance of  ln(S/R)', 'FontSize', 12);
    title({sprintf('[%s] Initial variogram parameter estimation', season_label), ...
           sprintf('[%s]  Starting values refined by WLS', mode_str)}, ...
        'FontWeight', 'bold', 'FontSize', 11);
    legend('Location', 'northeast', 'FontSize', 8.5);
    ylim([0, y_range * 1.4]);
    xlim([0, h_max/1000 * 1.02]);
    grid on; box on;

    saveas(fig0, fullfile(out_dir, 'Fig_B0_InitialParameters.png'));
    savefig(fig0, fullfile(out_dir, 'Fig_B0_InitialParameters.fig'));

    %% =================================================================
    %% PART 8 -- WLS fit, all 3 models (identical)
    %% =================================================================
    p_fit = cell(3,1); SSE = nan(3,1); WSSE = nan(3,1);
    w = sqrt(pair_n);

    fprintf('\n--- WLS fitting results (%s) ---\n', season_label);
    for m = 1:3
        res_fn = @(p) w .* (model_list{m}(p, pair_h) - pair_g);
        try
            p_fit{m} = lsqnonlin(res_fn, p0, lb, ub, opts);
        catch ME
            fprintf('Warning: WLS failed for %s (%s), using initial guess\n', model_names{m}, ME.message);
            p_fit{m} = p0;
        end
        resid   = model_list{m}(p_fit{m}, pair_h) - pair_g;
        SSE(m)  = sum(resid.^2);
        WSSE(m) = sum((w .* resid).^2);
        fprintf('%s:  C0=%.4f  CH=%.4f  AH=%6.1f km  SSE=%.4f  WSSE=%.4f\n', ...
            model_names{m}, p_fit{m}(1), p_fit{m}(2), p_fit{m}(3)/1000, SSE(m), WSSE(m));
    end
    [~, best_idx] = min(WSSE);
    fprintf('>>> Best model (%s): %s  (WSSE=%.4f)\n', season_label, model_names{best_idx}, WSSE(best_idx));

    %% =================================================================
    %% PART 9 -- Fig_B1: all 3 models on one plot (identical)
    %% =================================================================
    styles3 = {'-', '--', '-.'};
    fig1 = figure('Position', [50 50 860 560]);
    hold on;
    sc1 = scatter(plot_h/1000, plot_g, sz_display, 'k', 'filled', 'DisplayName', scatter_label);
    sc1.UserData = plot_n;
    if ~USE_BINS, sc1.MarkerFaceAlpha = 0.4; sc1.MarkerEdgeAlpha = 0; end

    for m = 1:3
        lw = 2.2 + (m == best_idx)*1.0;
        plot(h_plot/1000, model_list{m}(p_fit{m}, h_plot), styles3{m}, 'Color', colors3{m}, 'LineWidth', lw, ...
            'DisplayName', sprintf('%s   C_0=%.3f  C_H=%.3f  A_H=%.0f km   WSSE=%.4f%s', ...
                model_names{m}, p_fit{m}(1), p_fit{m}(2), p_fit{m}(3)/1000, WSSE(m), ...
                repmat('  *',1, m==best_idx)));
    end

    xlabel('Lag distance  h  (km)', 'FontSize', 12);
    ylabel('\gamma(h)  --  variance of  ln(S/R)', 'FontSize', 12);
    title({sprintf('[%s] Variogram model selection -- best fit: %s  *  [%s]', ...
            season_label, model_names{best_idx}, mode_str), ...
        'Parameters by Weighted Least Squares  (Cressie, 1985),  weight = sqrt(n_{pairs})'}, ...
        'FontWeight', 'bold', 'FontSize', 11);
    legend('Location', 'southeast', 'FontSize', 9);
    ylim([0, y_range * 1.45]);
    xlim([0, h_max/1000 * 1.02]);
    grid on; box on;

    saveas(fig1, fullfile(out_dir, 'Fig_B1_AllModels.png'));
    savefig(fig1, fullfile(out_dir, 'Fig_B1_AllModels.fig'));

    %% =================================================================
    %% PART 10 -- Fig_B2: three separate panels (identical)
    %% =================================================================
    fig2 = figure('Position', [50 50 1400 520]);
    for m = 1:3
        subplot(1,3,m); hold on;
        sc2 = scatter(plot_h/1000, plot_g, sz_display*0.7, 'k', 'filled', 'DisplayName', 'Empirical');
        sc2.UserData = plot_n;
        if ~USE_BINS, sc2.MarkerFaceAlpha = 0.4; sc2.MarkerEdgeAlpha = 0; end

        if USE_BINS
            for k = 1:numel(plot_h)
                g_fit = model_list{m}(p_fit{m}, plot_h(k));
                plot([plot_h(k)/1000 plot_h(k)/1000], [plot_g(k) g_fit], '-', ...
                    'Color', [0.75 0.75 0.75], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end

        plot(h_plot/1000, model_list{m}(p_fit{m}, h_plot), '-', 'Color', colors3{m}, 'LineWidth', 2.8, ...
            'DisplayName', model_names{m});
        plot(h_plot/1000, model_list{m}(p0, h_plot), ':', 'Color', colors3{m}*0.5+0.5, 'LineWidth', 1.2, ...
            'DisplayName', 'Initial guess');

        if m == best_idx
            title(sprintf('%s\nC_0=%.4f  C_H=%.4f  A_H=%.0f km\nWSSE=%.4f   * BEST', ...
                model_names{m}, p_fit{m}(1), p_fit{m}(2), p_fit{m}(3)/1000, WSSE(m)), ...
                'FontWeight','bold','FontSize',10);
        else
            title(sprintf('%s\nC_0=%.4f  C_H=%.4f  A_H=%.0f km\nWSSE=%.4f', ...
                model_names{m}, p_fit{m}(1), p_fit{m}(2), p_fit{m}(3)/1000, WSSE(m)), 'FontSize',10);
        end
        xlabel('h  (km)'); ylabel('\gamma(h)');
        ylim([0, y_range * 1.4]); xlim([0, h_max/1000 * 1.02]);
        legend('Location','southeast','FontSize',8);
        grid on; box on;
    end

    if USE_BINS
        sticks_note = 'Grey sticks = residuals (data minus fit)';
    else
        sticks_note = sprintf('All %d station pairs shown as individual points', n_pairs);
    end
    sgtitle({sprintf('[%s] Variogram model comparison', season_label), sticks_note, ...
        'Dotted line = initial guess,  Solid line = WLS optimized'}, 'FontWeight','bold');

    saveas(fig2, fullfile(out_dir, 'Fig_B2_ThreePanels.png'));
    savefig(fig2, fullfile(out_dir, 'Fig_B2_ThreePanels.fig'));

    %% =================================================================
    %% PART 11 -- Fig_B3: WSSE comparison bar chart (identical)
    %% =================================================================
    fig3 = figure('Position', [50 50 560 440]);
    bh_bar = bar(1:3, WSSE, 0.55);
    bh_bar.FaceColor = 'flat';
    for m = 1:3, bh_bar.CData(m,:) = colors3{m}; end
    hold on;
    bh2 = bar(best_idx, WSSE(best_idx), 0.55);
    bh2.FaceColor = colors3{best_idx};
    bh2.EdgeColor = 'k';
    bh2.LineWidth = 2.5;

    set(gca,'XTick',1:3,'XTickLabel',model_names,'FontSize',11);
    ylabel('Weighted SSE  (Cressie, 1985)','FontSize',11);
    title({sprintf('[%s] Variogram model fit quality', season_label), ...
        'Lower WSSE = better fit to empirical variogram'}, 'FontWeight','bold');
    ylim([0, max(WSSE)*1.18]);
    grid on; box on;

    saveas(fig3, fullfile(out_dir, 'Fig_B3_WSSE_Comparison.png'));
    savefig(fig3, fullfile(out_dir, 'Fig_B3_WSSE_Comparison.fig'));

    %% =================================================================
    %% PART 11b -- Fig_B3b: pair_n histogram (identical)
    %% =================================================================
    fig_hist = figure('Position', [50 50 700 460]);
    hold on;
    histogram(pair_n, 30, 'FaceColor', [0.20 0.47 0.75], 'EdgeColor', 'w', 'LineWidth', 0.8);
    xline(min(pair_n),    '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'DisplayName', sprintf('Min  = %d', min(pair_n)));
    xline(median(pair_n), '-',  'Color', [0.47 0.67 0.19], 'LineWidth', 2.0, 'DisplayName', sprintf('Median = %d', median(pair_n)));
    xline(mean(pair_n),   '--', 'Color', [0.93 0.69 0.13], 'LineWidth', 1.8, 'DisplayName', sprintf('Mean = %.0f', mean(pair_n)));
    xline(max(pair_n),    '--', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.8, 'DisplayName', sprintf('Max  = %d', max(pair_n)));

    xlabel('Number of shared valid hours  (n_{pairs})', 'FontSize', 12);
    ylabel('Count of station pairs', 'FontSize', 12);
    title({sprintf('[%s] Distribution of shared sample size across %d station pairs', season_label, n_pairs), ...
        sprintf('Min=%d   Median=%d   Mean=%.0f   Max=%d', min(pair_n), median(pair_n), mean(pair_n), max(pair_n))}, ...
        'FontWeight', 'bold', 'FontSize', 11);
    legend('Location', 'northeast', 'FontSize', 10);
    grid on; box on;

    saveas(fig_hist, fullfile(out_dir, 'Fig_B3b_PairSampleHistogram.png'));
    savefig(fig_hist, fullfile(out_dir, 'Fig_B3b_PairSampleHistogram.fig'));

    %% =================================================================
    %% PART 11c -- Fig_B3c: pair_n vs distance (identical)
    %% =================================================================
    fig_nh = figure('Position', [50 50 720 460]);
    hold on;
    scatter(pair_h/1000, pair_n, 18, pair_n, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0);
    colormap(gca, parula);
    cb = colorbar; cb.Label.String = 'n shared hours'; cb.Label.FontSize = 10;

    yline(min(pair_n),    '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', sprintf('Min  = %d h', min(pair_n)));
    yline(median(pair_n), '-',  'Color', [0.20 0.47 0.75], 'LineWidth', 1.8, 'DisplayName', sprintf('Median = %d h', median(pair_n)));
    yline(max(pair_n),    '--', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5, 'DisplayName', sprintf('Max  = %d h', max(pair_n)));

    xlabel('Distance between station pair  (km)', 'FontSize', 12);
    ylabel('Shared valid hours  (StRain>0  &  RaRain>0)', 'FontSize', 12);
    title({sprintf('[%s] Shared sample size per station pair  (%d pairs)', season_label, n_pairs), ...
        sprintf('Min=%d   Median=%d   Mean=%.0f   Max=%d  hours', min(pair_n), median(pair_n), mean(pair_n), max(pair_n))}, ...
        'FontWeight', 'bold', 'FontSize', 11);
    legend('Location', 'northeast', 'FontSize', 10);
    grid on; box on;

    saveas(fig_nh, fullfile(out_dir, 'Fig_B3c_PairN_vs_Distance.png'));
    savefig(fig_nh, fullfile(out_dir, 'Fig_B3c_PairN_vs_Distance.fig'));

    %% =================================================================
    %% PART 12 -- Fig_B4: final selected model (identical)
    %% =================================================================
    fig4 = figure('Position', [50 50 720 520]);
    hold on;
    sc4 = scatter(plot_h/1000, plot_g, sz_display, 'k', 'filled', 'DisplayName', scatter_label);
    sc4.UserData = plot_n;
    if ~USE_BINS, sc4.MarkerFaceAlpha = 0.4; sc4.MarkerEdgeAlpha = 0; end

    plot(h_plot/1000, model_list{best_idx}(p_fit{best_idx}, h_plot), '-', ...
        'Color', colors3{best_idx}, 'LineWidth', 3.2, 'DisplayName', sprintf('%s model', model_names{best_idx}));

    yline(p_fit{best_idx}(1), '--', 'Color',[0.85 0.33 0.10], 'LineWidth',1.5, ...
        'DisplayName', sprintf('Nugget  C_0 = %.4f', p_fit{best_idx}(1)));
    yline(p_fit{best_idx}(1)+p_fit{best_idx}(2), '--', 'Color',[0.47 0.67 0.19], 'LineWidth',1.5, ...
        'DisplayName', sprintf('Sill  C_0+C_H = %.4f', p_fit{best_idx}(1)+p_fit{best_idx}(2)));
    xline(p_fit{best_idx}(3)/1000, '--', 'Color',[0.20 0.47 0.75], 'LineWidth',1.5, ...
        'DisplayName', sprintf('Range  A_H = %.0f km', p_fit{best_idx}(3)/1000));

    xlabel('Lag distance  h  (km)', 'FontSize', 12);
    ylabel('\gamma(h)  --  variance of  ln(S/R)', 'FontSize', 12);
    title({sprintf('[%s] Selected variogram model: %s', season_label, model_names{best_idx}), ...
        sprintf('C_0=%.4f   C_H=%.4f   A_H=%.0f km   WSSE=%.4f', ...
            p_fit{best_idx}(1), p_fit{best_idx}(2), p_fit{best_idx}(3)/1000, WSSE(best_idx))}, ...
        'FontWeight','bold','FontSize',11);
    legend('Location','southeast','FontSize',10);
    ylim([0, y_range * 1.35]);
    xlim([0, h_max/1000 * 1.02]);
    grid on; box on;

    saveas(fig4, fullfile(out_dir, 'Fig_B4_SelectedModel_Final.png'));
    savefig(fig4, fullfile(out_dir, 'Fig_B4_SelectedModel_Final.fig'));

    %% =================================================================
    %% PART 12b -- custom tooltips (identical)
    %% =================================================================
    cursor_fn = @(~, event) custom_variogram_tooltip(event, plot_h, plot_g, plot_n);
    for fig_h = [fig0, fig1, fig4]
        dcm = datacursormode(fig_h);
        dcm.UpdateFcn = cursor_fn;
        dcm.Enable = 'on';
    end

    %% =================================================================
    %% PART 13 -- tables (identical)
    %% =================================================================
    T_fit = table(model_names', ...
        [p0(1); p0(1); p0(1)], [p0(2); p0(2); p0(2)], [p0(3)/1000; p0(3)/1000; p0(3)/1000], ...
        cellfun(@(p) p(1), p_fit), cellfun(@(p) p(2), p_fit), cellfun(@(p) p(3)/1000, p_fit), ...
        SSE, WSSE, ...
        'VariableNames',{'Model','Init_C0','Init_CH','Init_AH_km', ...
            'WLS_C0','WLS_CH','WLS_AH_km','SSE','WSSE'});
    writetable(T_fit, fullfile(out_dir,'Table_VariogramFit.csv'));

    T_init = table( ...
        {'Nugget C0';'Partial sill CH';'Total sill';'Range AH';'Source threshold';'Display mode';'Season'}, ...
        {sprintf('%.4f', C0_init); sprintf('%.4f', CH_init); sprintf('%.4f', C0_init+CH_init); ...
         sprintf('%.1f km', AH_init/1000); sprintf('h < %.1f km (10th pct)', short_threshold/1000); ...
         mode_str; season_label}, ...
        {'Mean gamma for h below the 10th-pct distance threshold'; '85th pct of gamma minus C0'; ...
         'C0 + CH'; 'h where gamma = C0 + 0.632*CH'; '10th percentile of all inter-station distances'; ...
         'Controlled by USE_BINS flag'; 'Season this table corresponds to'}, ...
        'VariableNames',{'Parameter','Value','Method'});
    writetable(T_init, fullfile(out_dir,'Table_InitialEstimates.csv'));

    fprintf('\n--- Final summary (%s) ---\n', season_label);
    fprintf('Best model:   %s\n', model_names{best_idx});
    fprintf('  C0=%.4f  CH=%.4f  AH=%.1f km  WSSE=%.4f\n', ...
        p_fit{best_idx}(1), p_fit{best_idx}(2), p_fit{best_idx}(3)/1000, WSSE(best_idx));
    fprintf('Saved to: %s\n', out_dir);

    close([fig0, fig1, fig2, fig3, fig_hist, fig_nh, fig4]);   % keep memory clean across the 4 iterations
end

fprintf('\n\nAll four seasons complete. Each has its own output folder with the full 7-figure, 2-table pipeline.\n');

%% =====================================================================
%% helper function (identical to the uploaded script)
%% =====================================================================

%% ============================================================
%  Local functions (moved to the end of the script so it also runs
%  on MATLAB releases older than R2024a)
%% ============================================================
function txt = custom_variogram_tooltip(event, ph, pg, pn)
pos = event.Position;
dx = (ph/1000 - pos(1)) ./ range(ph/1000);
dy = (pg      - pos(2)) ./ range(pg);
[~, idx] = min(dx.^2 + dy.^2);
txt = {
    sprintf('h  = %.2f km',  ph(idx)/1000), ...
    sprintf('gamma  = %.5f', pg(idx)),       ...
    sprintf('n shared hours = %d', pn(idx))  ...
};
end
