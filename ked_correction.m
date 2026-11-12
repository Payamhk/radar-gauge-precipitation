%% ------------------------------------------------------------------
%  ked_correction.m
%  Part of: radar-gauge-precipitation (Payam Heidarian, University of
%  Brescia / MeteoSwiss collaboration). Inline comments are partly in
%  Persian. Data are not included; set the <DATA_ROOT> paths below to
%  your own folders before running.
%% ------------------------------------------------------------------
%% ============================================================
%  VALIDATION + CORRECTION OUTPUT (KED VERSION) -- CSV6 ایستگاه‌ها
%  ------------------------------------------------------------
%  این نسخه دقیقاً همان خروجی‌های اسکریپت OK را تولید می‌کند؛ تنها
%  تفاوت، موتور کریجینگ است:
%    OK : هدف = ln(S/R)   ->  Corrected = Ra * exp(pred)
%    KED: هدف = ln(S)،  درفت = ln(R)   ->  Corrected = exp(pred)
%  پارامترهای واریوگرام (فصلی + سالانه) عیناً از نسخهٔ OK قفل شده‌اند.
%
%  *** این نسخه: Safeguard 1 و 2 قابل‌تنظیم، دقیقاً مثل نسخهٔ ARPA ***
%
%  منطق فیلتر/fallback/معیارها هیچ تغییری نکرده:
%  - همه ردیف‌های هر CSV6 خوانده می‌شوند (بدون فیلتر اولیه)
%  - کریجینگ: فقط ساعت‌هایی که StRain > 0 و RaRain > 0
%  - در غیر این صورت: همون RaRain خام (fallback)
%  - NRMSE/r/RD%: فقط ساعت‌هایی که RaRain > 0 (CPC_orig > 0)
%  - خروجی CSV: همه ستون‌های اصلی + CorrectedRain_Seasonal_KED + CorrectedRain_Annual_KED
%  - گرید زمانی: اتحاد بازه train و CSV6 (extended grid)
%% ============================================================
clc; clear;

%% ============================================================
%  تنظیمات Safeguard -- دلخواه شما (همان الگوی نسخهٔ ARPA)
%  ------------------------------------------------------------
%  USE_SAFEGUARD_1: رد کردن ساعت‌هایی که target_lnR بیرون از بازهٔ
%                   ln(R) همسایه‌هاست (برون‌یابی در درفت)
%  USE_SAFEGUARD_2: رد کردن ساعت‌هایی که max(|lambda|) از آستانهٔ
%                   زیر بیشتر شود
%  SAFEGUARD_2_THRESHOLD: آستانهٔ Safeguard 2 (فقط وقتی فعال است
%                   اثر دارد)
%% ============================================================
USE_SAFEGUARD_1       = true;
USE_SAFEGUARD_2       = true;
SAFEGUARD_2_THRESHOLD = 10;

path_gs = "<DATA_ROOT>\" + ...
    "Profession\1-PRECIPITATION\Code\CODE\Moghayese gs vs nc franciacorta\" + ...
    "چیزی که استاد خواسته\code\1-Combi-GS-Valid H\1-in\GS-cel-xy-val H-pre";

cfg_tag = sprintf('SG1_%d_SG2_%d_Thr_%g', USE_SAFEGUARD_1, USE_SAFEGUARD_2, SAFEGUARD_2_THRESHOLD);
out_dir = ['Output_CSV6_Corrected_KED_' cfg_tag];
if ~exist(out_dir,'dir'), mkdir(out_dir); end

csv6_names = {
    'Consorzio_Clusane_Belvedere_6620_corretto'
    'Consorzio_Sant_Eusebio_40485_corretto'
    'Consorzio_Brescia_Sant_Anna_6619_corretto'
    'Consorzio_Dosso_Baiana_40484_corretto'
    'Consorzio_Crociate_40479_corretto'
    'Consorzio_Versella_40480_corretto'
};
n6 = numel(csv6_names);

season_map   = [4,4,1,1,1,2,2,2,3,3,3,4];
month_names  = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

%% ============================================================
%  پارامترهای واریوگرام (locked) -- عیناً از نسخهٔ OK
%  radius تابستان(2) و پاییز(3) = 10000 | بهار(1) و زمستان(4) = Inf
%% ============================================================
vario_seasonal(1).C0=0.1008; vario_seasonal(1).CH=0.7552; vario_seasonal(1).AH=19000; vario_seasonal(1).radius=Inf;
vario_seasonal(2).C0=0.0743; vario_seasonal(2).CH=0.4099; vario_seasonal(2).AH= 8000; vario_seasonal(2).radius=10000;
vario_seasonal(3).C0=0.0286; vario_seasonal(3).CH=0.4147; vario_seasonal(3).AH= 9000; vario_seasonal(3).radius=10000;
vario_seasonal(4).C0=0.0258; vario_seasonal(4).CH=0.5605; vario_seasonal(4).AH=18000; vario_seasonal(4).radius=Inf;

vario_annual.C0     = 0.0561;
vario_annual.CH     = 0.5016;
vario_annual.AH     = 12000;
vario_annual.radius = Inf;

exp_model    = @(C0,CH,AH,h) C0 + CH.*(1 - exp(-h./AH));
min_stations = 3;

%% ============================================================
%  خواندن ایستگاه‌های train -- اینجا lnS و lnR جدا نگه داشته می‌شوند
%% ============================================================
all_files   = dir(fullfile(path_gs,'Consorzio*.csv'));
all_names   = erase({all_files.name},'.csv');
train_names = setdiff(all_names, csv6_names,'stable');

raw_train = cell(numel(train_names),1);
cnt = 0;
for i = 1:numel(train_names)
    fpath = fullfile(path_gs,[train_names{i} '.csv']);
    if ~isfile(fpath), continue; end
    T = readtable(fpath);
    T = T(:,1:7);
    T.Properties.VariableNames = {'Time','Sx','Sy','Rx','Ry','StRain','RaRain'};
    T.Time = datetime(T.Time,'InputFormat','yyyy-MM-dd HH:mm:ss');
    valid  = T.StRain > 0 & T.RaRain > 0;
    if sum(valid) < 2, continue; end
    cnt = cnt+1;
    raw_train{cnt}.Sx   = T.Sx(1);
    raw_train{cnt}.Sy   = T.Sy(1);
    raw_train{cnt}.Time = T.Time(valid);
    raw_train{cnt}.lnS  = log(T.StRain(valid));
    raw_train{cnt}.lnR  = log(T.RaRain(valid));
end
raw_train = raw_train(1:cnt);
n_train   = cnt;

Sx_tr = cellfun(@(r) r.Sx, raw_train)';
Sy_tr = cellfun(@(r) r.Sy, raw_train)';

train_start = min(cellfun(@(r) min(r.Time), raw_train));
train_end   = max(cellfun(@(r) max(r.Time), raw_train));
fprintf('Training range: %s  to  %s\n', train_start, train_end);

%% ============================================================
%  پاس سریع روی CSV6 برای بازه زمانی کامل
%% ============================================================
csv6_min_dt = datetime.empty(0,1);
csv6_max_dt = datetime.empty(0,1);
for v = 1:n6
    fpath = fullfile(path_gs,[csv6_names{v} '.csv']);
    if ~isfile(fpath), continue; end
    Ttmp = readtable(fpath);
    Ttmp.Properties.VariableNames(1) = {'Time'};
    Ttmp.Time = datetime(Ttmp.Time,'InputFormat','yyyy-MM-dd HH:mm:ss');
    csv6_min_dt(end+1,1) = min(Ttmp.Time); %#ok<AGROW>
    csv6_max_dt(end+1,1) = max(Ttmp.Time); %#ok<AGROW>
end
csv6_start = min(csv6_min_dt);
csv6_end   = max(csv6_max_dt);
fprintf('CSV6 full range: %s  to  %s\n', csv6_start, csv6_end);

%% ============================================================
%  گرید زمانی گسترش‌یافته
%% ============================================================
grid_start = min(train_start, csv6_start);
grid_end   = max(train_end,   csv6_end);
grid_times = (grid_start:hours(1):grid_end)';
nGrid      = numel(grid_times);
season_of_grid = season_map(month(grid_times))';
month_of_grid  = month(grid_times);
fprintf('Extended grid: %s  to  %s  (%d hours)\n', grid_start, grid_end, nGrid);

lnS_dense_tr = nan(nGrid, n_train);
lnR_dense_tr = nan(nGrid, n_train);
for s = 1:n_train
    idx = round(hours(raw_train{s}.Time - grid_start)) + 1;
    ok  = idx >= 1 & idx <= nGrid;
    lnS_dense_tr(idx(ok), s) = raw_train{s}.lnS(ok);
    lnR_dense_tr(idx(ok), s) = raw_train{s}.lnR(ok);
end
Dfull_tr = squareform(pdist([Sx_tr(:), Sy_tr(:)]));

%% ============================================================
%  کریجینگ روی همه CSV6 (با موتور KED، Safeguard قابل‌تنظیم)
%% ============================================================
St_mat    = nan(nGrid, n6);
Ra_mat    = nan(nGrid, n6);
Corr_seas = nan(nGrid, n6);
Corr_ann  = nan(nGrid, n6);

csv6_data = cell(n6,1);

diag_extrap_seas = []; diag_maxlam_seas = [];
diag_extrap_ann  = []; diag_maxlam_ann  = [];

for v = 1:n6
    fpath = fullfile(path_gs,[csv6_names{v} '.csv']);
    if ~isfile(fpath)
        fprintf('File not found: %s\n', csv6_names{v});
        continue;
    end

    T_raw = readtable(fpath);
    T_raw = T_raw(:,1:7);
    T_raw.Properties.VariableNames = {'Time','Sx','Sy','Rx','Ry','StRain','RaRain'};
    T_raw.Time = datetime(T_raw.Time,'InputFormat','yyyy-MM-dd HH:mm:ss');

    nRows = height(T_raw);
    Sx_v  = T_raw.Sx(1);
    Sy_v  = T_raw.Sy(1);
    fprintf('Processing %s  (%d rows) ...\n', csv6_names{v}, nRows);

    CorrectedRain_Seasonal_KED = T_raw.RaRain;
    CorrectedRain_Annual_KED   = T_raw.RaRain;

    for r = 1:nRows
        t_idx = round(hours(T_raw.Time(r) - grid_start)) + 1;
        if t_idx < 1 || t_idx > nGrid, continue; end

        if ~isnan(T_raw.StRain(r)) && ~isnan(T_raw.RaRain(r))
            St_mat(t_idx,v) = T_raw.StRain(r);
            Ra_mat(t_idx,v) = T_raw.RaRain(r);
        end

        if ~isnan(T_raw.RaRain(r))
            Corr_seas(t_idx,v) = T_raw.RaRain(r);
            Corr_ann(t_idx,v)  = T_raw.RaRain(r);
        end

        if isnan(T_raw.StRain(r)) || isnan(T_raw.RaRain(r)) || ...
           T_raw.StRain(r) <= 0   || T_raw.RaRain(r) <= 0
            continue;
        end

        s_idx = season_of_grid(t_idx);

        [pred_seas, fx_s, ml_s] = do_kriging_KED(t_idx, Sx_v, Sy_v, ...
            vario_seasonal(s_idx).C0, vario_seasonal(s_idx).CH, ...
            vario_seasonal(s_idx).AH, vario_seasonal(s_idx).radius, ...
            lnS_dense_tr, lnR_dense_tr, Dfull_tr, Sx_tr, Sy_tr, n_train, ...
            T_raw.RaRain(r), exp_model, min_stations, ...
            USE_SAFEGUARD_1, USE_SAFEGUARD_2, SAFEGUARD_2_THRESHOLD);

        [pred_ann, fx_a, ml_a] = do_kriging_KED(t_idx, Sx_v, Sy_v, ...
            vario_annual.C0, vario_annual.CH, ...
            vario_annual.AH, vario_annual.radius, ...
            lnS_dense_tr, lnR_dense_tr, Dfull_tr, Sx_tr, Sy_tr, n_train, ...
            T_raw.RaRain(r), exp_model, min_stations, ...
            USE_SAFEGUARD_1, USE_SAFEGUARD_2, SAFEGUARD_2_THRESHOLD);

        diag_extrap_seas(end+1) = fx_s; %#ok<AGROW>
        diag_extrap_ann(end+1)  = fx_a; %#ok<AGROW>

        if ~isnan(pred_seas)
            diag_maxlam_seas(end+1) = ml_s; %#ok<AGROW>
            CorrectedRain_Seasonal_KED(r) = pred_seas;
            Corr_seas(t_idx,v)            = pred_seas;
        end

        if ~isnan(pred_ann)
            diag_maxlam_ann(end+1) = ml_a; %#ok<AGROW>
            CorrectedRain_Annual_KED(r) = pred_ann;
            Corr_ann(t_idx,v)           = pred_ann;
        end
    end

    T_out = T_raw;
    T_out.CorrectedRain_Seasonal_KED = CorrectedRain_Seasonal_KED;
    T_out.CorrectedRain_Annual_KED   = CorrectedRain_Annual_KED;
    csv6_data{v} = T_out;

    safe_name = matlab.lang.makeValidName(csv6_names{v});
    out_fname  = fullfile(out_dir, sprintf('%s_corrected_KED.csv', safe_name));
    writetable(T_out, out_fname);
    fprintf('  Saved CSV: %s\n', out_fname);
end

fprintf('\nKED kriging done. [Safeguard 1 = %s | Safeguard 2 = %s (threshold=%g)]\n', ...
    mat2str(USE_SAFEGUARD_1), mat2str(USE_SAFEGUARD_2), SAFEGUARD_2_THRESHOLD);

%% ============================================================
%  خلاصهٔ تشخیصی (CSV6)
%% ============================================================
fprintf('\n=== DIAGNOSTIC SUMMARY (CSV6) ===\n');
fprintf('Config: Safeguard1=%s | Safeguard2=%s | Threshold2=%g\n', ...
    mat2str(USE_SAFEGUARD_1), mat2str(USE_SAFEGUARD_2), SAFEGUARD_2_THRESHOLD);
fprintf('--- Seasonal model ---\n');
fprintf('  N attempted (total hours kriging was tried): %d\n', numel(diag_extrap_seas));
fprintf('  N accepted (passed active safeguards): %d\n', numel(diag_maxlam_seas));
fprintf('  Extrapolation rate (would trigger Safeguard 1 if enabled): %.1f%%\n', 100*mean(diag_extrap_seas));
fprintf('  max(|lambda|) percentiles -> p50=%.2f  p90=%.2f  p95=%.2f  p99=%.2f  p99.9=%.2f  max=%.2f\n', ...
    prctile(diag_maxlam_seas,50), prctile(diag_maxlam_seas,90), ...
    prctile(diag_maxlam_seas,95), prctile(diag_maxlam_seas,99), ...
    prctile(diag_maxlam_seas,99.9), max(diag_maxlam_seas));
fprintf('--- Annual model ---\n');
fprintf('  N attempted (total hours kriging was tried): %d\n', numel(diag_extrap_ann));
fprintf('  N accepted (passed active safeguards): %d\n', numel(diag_maxlam_ann));
fprintf('  Extrapolation rate (would trigger Safeguard 1 if enabled): %.1f%%\n', 100*mean(diag_extrap_ann));
fprintf('  max(|lambda|) percentiles -> p50=%.2f  p90=%.2f  p95=%.2f  p99=%.2f  p99.9=%.2f  max=%.2f\n', ...
    prctile(diag_maxlam_ann,50), prctile(diag_maxlam_ann,90), ...
    prctile(diag_maxlam_ann,95), prctile(diag_maxlam_ann,99), ...
    prctile(diag_maxlam_ann,99.9), max(diag_maxlam_ann));

fig_diag = figure('Position',[50 50 1000 420]);
subplot(1,2,1);
histogram(diag_maxlam_seas, 50, 'FaceColor',[0.20 0.47 0.75]);
xlabel('max(|\lambda|)  per hour'); ylabel('count');
title('Seasonal model -- weight magnitude distribution');
grid on; box on;
subplot(1,2,2);
histogram(diag_maxlam_ann, 50, 'FaceColor',[0.47 0.67 0.19]);
xlabel('max(|\lambda|)  per hour'); ylabel('count');
title('Annual model -- weight magnitude distribution');
grid on; box on;
sgtitle('CSV6 -- distribution of KED weight magnitudes','FontWeight','bold');
saveas(fig_diag, fullfile(out_dir,'Fig_Diag_LambdaDistribution.png'));

%% ============================================================
%  NRMSE ماهانه -- pooled  (فیلتر: RaRain > 0)
%% ============================================================
NRMSE_ra   = nan(12,1);
NRMSE_seas = nan(12,1);
NRMSE_ann  = nan(12,1);
N_month    = zeros(12,1);

for m = 1:12
    t_mask = month_of_grid == m;
    St_p = St_mat(t_mask,:);
    Ra_p = Ra_mat(t_mask,:);
    Cs_p = Corr_seas(t_mask,:);
    Ca_p = Corr_ann(t_mask,:);

    ok = ~isnan(St_p) & ~isnan(Ra_p) & Ra_p > 0 & ~isnan(Cs_p) & ~isnan(Ca_p);

    N_month(m) = sum(ok(:));
    if N_month(m) < 5, continue; end

    St_v = St_p(ok); Ra_v = Ra_p(ok); Cs_v = Cs_p(ok); Ca_v = Ca_p(ok);
    nv   = mean(St_v);
    NRMSE_ra(m)   = sqrt(mean((Ra_v - St_v).^2)) / nv;
    NRMSE_seas(m) = sqrt(mean((Cs_v - St_v).^2)) / nv;
    NRMSE_ann(m)  = sqrt(mean((Ca_v - St_v).^2)) / nv;
end

ok_a  = ~isnan(St_mat) & ~isnan(Ra_mat) & Ra_mat > 0 & ~isnan(Corr_seas) & ~isnan(Corr_ann);
St_a  = St_mat(ok_a); Ra_a = Ra_mat(ok_a);
Cs_a  = Corr_seas(ok_a); Ca_a = Corr_ann(ok_a);
nv_a  = mean(St_a);
NRMSE_ann_ra   = sqrt(mean((Ra_a - St_a).^2)) / nv_a;
NRMSE_ann_seas = sqrt(mean((Cs_a - St_a).^2)) / nv_a;
NRMSE_ann_ann  = sqrt(mean((Ca_a - St_a).^2)) / nv_a;

fprintf('\n--- Annual NRMSE (CSV6, KED, pooled, filter: RaRain>0) ---\n');
fprintf('Radar:    %.4f\n', NRMSE_ann_ra);
fprintf('Seasonal: %.4f  (Delta=%.1f%%)\n', NRMSE_ann_seas, 100*(NRMSE_ann_ra-NRMSE_ann_seas)/NRMSE_ann_ra);
fprintf('Annual:   %.4f  (Delta=%.1f%%)\n', NRMSE_ann_ann,  100*(NRMSE_ann_ra-NRMSE_ann_ann)/NRMSE_ann_ra);

%% ============================================================
%  Pearson r ماهانه -- pooled
%% ============================================================
R_ra   = nan(12,1);
R_seas = nan(12,1);
R_ann  = nan(12,1);

for m = 1:12
    t_mask = month_of_grid == m;
    St_p = St_mat(t_mask,:);
    Ra_p = Ra_mat(t_mask,:);
    Cs_p = Corr_seas(t_mask,:);
    Ca_p = Corr_ann(t_mask,:);

    ok = ~isnan(St_p) & ~isnan(Ra_p) & Ra_p > 0 & ~isnan(Cs_p) & ~isnan(Ca_p);
    if sum(ok(:)) < 5, continue; end

    St_v = St_p(ok); Ra_v = Ra_p(ok); Cs_v = Cs_p(ok); Ca_v = Ca_p(ok);
    tmp = corrcoef(Ra_v, St_v);   R_ra(m)   = tmp(1,2);
    tmp = corrcoef(Cs_v, St_v);   R_seas(m) = tmp(1,2);
    tmp = corrcoef(Ca_v, St_v);   R_ann(m)  = tmp(1,2);
end

tmp = corrcoef(Ra_a, St_a); R_ann_ra   = tmp(1,2);
tmp = corrcoef(Cs_a, St_a); R_ann_seas = tmp(1,2);
tmp = corrcoef(Ca_a, St_a); R_ann_ann  = tmp(1,2);

fprintf('\n--- Annual Pearson r (CSV6, KED, pooled, filter: RaRain>0) ---\n');
fprintf('Radar:    %.4f\n', R_ann_ra);
fprintf('Seasonal: %.4f\n', R_ann_seas);
fprintf('Annual:   %.4f\n', R_ann_ann);

%% ============================================================
%  RD% ماهانه و سالانه -- pooled  (میانگین فضایی اول، بعد جمع)
%% ============================================================
RD_ra   = nan(12,1);
RD_seas = nan(12,1);
RD_ann  = nan(12,1);

sp_St_rd   = nan(nGrid,1);
sp_Ra_rd   = nan(nGrid,1);
sp_Seas_rd = nan(nGrid,1);
sp_Ann_rd  = nan(nGrid,1);

for t = 1:nGrid
    ok_t = ~isnan(St_mat(t,:)) & ~isnan(Ra_mat(t,:)) & Ra_mat(t,:) > 0 & ...
           ~isnan(Corr_seas(t,:)) & ~isnan(Corr_ann(t,:));
    if sum(ok_t) == 0, continue; end
    sp_St_rd(t)   = mean(St_mat(t,ok_t));
    sp_Ra_rd(t)   = mean(Ra_mat(t,ok_t));
    sp_Seas_rd(t) = mean(Corr_seas(t,ok_t));
    sp_Ann_rd(t)  = mean(Corr_ann(t,ok_t));
end

for m = 1:12
    t_mask = month_of_grid == m;
    St_m   = sp_St_rd(t_mask);
    Ra_m   = sp_Ra_rd(t_mask);
    Cs_m   = sp_Seas_rd(t_mask);
    Ca_m   = sp_Ann_rd(t_mask);

    ok_m = ~isnan(St_m) & ~isnan(Ra_m) & ~isnan(Cs_m) & ~isnan(Ca_m);
    if sum(ok_m) < 5, continue; end

    sum_gs     = sum(St_m(ok_m));
    RD_ra(m)   = (sum(Ra_m(ok_m))  - sum_gs) / sum_gs * 100;
    RD_seas(m) = (sum(Cs_m(ok_m))  - sum_gs) / sum_gs * 100;
    RD_ann(m)  = (sum(Ca_m(ok_m))  - sum_gs) / sum_gs * 100;
end

ok_rd_all   = ~isnan(sp_St_rd) & ~isnan(sp_Ra_rd) & ~isnan(sp_Seas_rd) & ~isnan(sp_Ann_rd);
sum_gs_all  = sum(sp_St_rd(ok_rd_all));
RD_ann_ra   = (sum(sp_Ra_rd(ok_rd_all))   - sum_gs_all) / sum_gs_all * 100;
RD_ann_seas = (sum(sp_Seas_rd(ok_rd_all)) - sum_gs_all) / sum_gs_all * 100;
RD_ann_ann  = (sum(sp_Ann_rd(ok_rd_all))  - sum_gs_all) / sum_gs_all * 100;

fprintf('\n--- Annual RD%% (CSV6, KED, spatial mean first then sum, filter: RaRain>0) ---\n');
fprintf('Radar:    %.2f%%\n', RD_ann_ra);
fprintf('Seasonal: %.2f%%\n', RD_ann_seas);
fprintf('Annual:   %.2f%%\n', RD_ann_ann);

%% ============================================================
%  کومولاتیو pooled
%% ============================================================
all_T_cum    = datetime.empty;
all_gs_cum   = [];
all_ra_cum   = [];
all_seas_cum = [];
all_ann_cum  = [];

for v = 1:n6
    for t = 1:nGrid
        if isnan(St_mat(t,v)) || isnan(Ra_mat(t,v)) || ...
           isnan(Corr_seas(t,v)) || isnan(Corr_ann(t,v))
            continue;
        end
        all_T_cum    = [all_T_cum;    grid_times(t)]; %#ok<AGROW>
        all_gs_cum   = [all_gs_cum;   St_mat(t,v)]; %#ok<AGROW>
        all_ra_cum   = [all_ra_cum;   Ra_mat(t,v)]; %#ok<AGROW>
        all_seas_cum = [all_seas_cum; Corr_seas(t,v)]; %#ok<AGROW>
        all_ann_cum  = [all_ann_cum;  Corr_ann(t,v)]; %#ok<AGROW>
    end
end

[uT_pool, ~, ic_pool] = unique(all_T_cum);
sp_gs_h   = accumarray(ic_pool, all_gs_cum,   [], @(x)mean(x,'omitnan'));
sp_ra_h   = accumarray(ic_pool, all_ra_cum,   [], @(x)mean(x,'omitnan'));
sp_seas_h = accumarray(ic_pool, all_seas_cum, [], @(x)mean(x,'omitnan'));
sp_ann_h  = accumarray(ic_pool, all_ann_cum,  [], @(x)mean(x,'omitnan'));

times_ok = uT_pool;
cumSt   = cumsum(sp_gs_h);
cumRa   = cumsum(sp_ra_h);
cumSeas = cumsum(sp_seas_h);
cumAnn  = cumsum(sp_ann_h);

relDiff_Ra   = (cumRa   - cumSt) ./ cumSt * 100;
relDiff_Seas = (cumSeas - cumSt) ./ cumSt * 100;
relDiff_Ann  = (cumAnn  - cumSt) ./ cumSt * 100;

fprintf('\nFinal cumulative (pooled, KED): Station=%.1f mm | Radar=%.1f mm | Seasonal=%.1f mm | Annual=%.1f mm\n', ...
    cumSt(end), cumRa(end), cumSeas(end), cumAnn(end));

%% ============================================================
%  نمودار ۱: NRMSE ماهانه pooled
%% ============================================================
valid_m = N_month >= 5;
x = find(valid_m);

fig1 = figure('Position',[50 50 1000 500]); hold on;
b1 = bar(x-0.25, NRMSE_ra(valid_m),   0.24, 'FaceColor',[0.6 0.6 0.6]);
b2 = bar(x,      NRMSE_seas(valid_m),  0.24, 'FaceColor',[0.20 0.47 0.75]);
b3 = bar(x+0.25, NRMSE_ann(valid_m),   0.24, 'FaceColor',[0.47 0.67 0.19]);
set(gca,'XTick',1:12,'XTickLabel',month_names,'FontSize',10);
ylabel('NRMSE','FontSize',11);
title({'Monthly NRMSE (KED) -- CSV6 (filter: RaRain>0)', ...
    sprintf('Annual -- Radar:%.3f  Seasonal KED:%.3f(%.0f%%)  Annual KED:%.3f(%.0f%%)', ...
        NRMSE_ann_ra, NRMSE_ann_seas, 100*(NRMSE_ann_ra-NRMSE_ann_seas)/NRMSE_ann_ra, ...
        NRMSE_ann_ann,  100*(NRMSE_ann_ra-NRMSE_ann_ann)/NRMSE_ann_ra)}, 'FontWeight','bold');
legend([b1 b2 b3],{'Radar (before)','Seasonal KED','Annual KED'},'Location','northwest','FontSize',9);
grid on; box on;
saveas(fig1, fullfile(out_dir,'Fig_V1_KED_NRMSE_Monthly_Comparison.png'));
savefig(fig1, fullfile(out_dir,'Fig_V1_KED_NRMSE_Monthly_Comparison.fig'));

%% ============================================================
%  نمودار ۱b: Pearson r ماهانه pooled
%% ============================================================
valid_r = ~isnan(R_ra) | ~isnan(R_seas) | ~isnan(R_ann);
xr = find(valid_r);

fig1b = figure('Position',[50 50 1000 500]); hold on;
b1r = bar(xr-0.25, R_ra(valid_r),   0.24, 'FaceColor',[0.6 0.6 0.6]);
b2r = bar(xr,      R_seas(valid_r),  0.24, 'FaceColor',[0.20 0.47 0.75]);
b3r = bar(xr+0.25, R_ann(valid_r),   0.24, 'FaceColor',[0.47 0.67 0.19]);
yline(1,'k--','LineWidth',1.2,'HandleVisibility','off');
set(gca,'XTick',1:12,'XTickLabel',month_names,'FontSize',10);
ylim([0 1.1]);
ylabel('Pearson r','FontSize',11);
title({'Monthly Pearson r (KED) -- CSV6 (filter: RaRain>0)', ...
    sprintf('Annual -- Radar:%.3f  Seasonal KED:%.3f  Annual KED:%.3f', ...
        R_ann_ra, R_ann_seas, R_ann_ann)}, 'FontWeight','bold');
legend([b1r b2r b3r],{'Radar (before)','Seasonal KED','Annual KED'},'Location','southwest','FontSize',9);
grid on; box on;
saveas(fig1b, fullfile(out_dir,'Fig_V1b_KED_PearsonR_Monthly_Comparison.png'));
savefig(fig1b, fullfile(out_dir,'Fig_V1b_KED_PearsonR_Monthly_Comparison.fig'));

%% ============================================================
%  نمودار ۱c: RD% ماهانه pooled
%% ============================================================
valid_rd = ~isnan(RD_ra) | ~isnan(RD_seas) | ~isnan(RD_ann);
xrd = find(valid_rd);

fig1c = figure('Position',[50 50 1000 500]); hold on;
b1c = bar(xrd-0.25, RD_ra(valid_rd),   0.24, 'FaceColor',[0.6 0.6 0.6]);
b2c = bar(xrd,      RD_seas(valid_rd),  0.24, 'FaceColor',[0.20 0.47 0.75]);
b3c = bar(xrd+0.25, RD_ann(valid_rd),   0.24, 'FaceColor',[0.47 0.67 0.19]);
yline(0,'k-','LineWidth',1.5,'HandleVisibility','off');
set(gca,'XTick',1:12,'XTickLabel',month_names,'FontSize',10);
ylabel('RD%  (CPC - GS) / GS x 100','FontSize',11);
title({'Monthly RD% (KED) -- CSV6 (filter: RaRain>0, cumulative sum per month)', ...
    sprintf('Annual -- Radar:%.1f%%  Seasonal KED:%.1f%%  Annual KED:%.1f%%', ...
        RD_ann_ra, RD_ann_seas, RD_ann_ann)}, 'FontWeight','bold');
legend([b1c b2c b3c],{'Radar (before)','Seasonal KED','Annual KED'},'Location','northeast','FontSize',9);
grid on; box on;
saveas(fig1c, fullfile(out_dir,'Fig_V1c_KED_RD_Monthly_Comparison.png'));
savefig(fig1c, fullfile(out_dir,'Fig_V1c_KED_RD_Monthly_Comparison.fig'));

%% ============================================================
%  نمودار ۱d: نتایج سالانه -- سه پانل کنار هم
%% ============================================================
methods     = {'Radar','Seasonal KED','Annual KED'};
clrs        = [0.6 0.6 0.6; 0.20 0.47 0.75; 0.47 0.67 0.19];
nrmse_vals  = [NRMSE_ann_ra;   NRMSE_ann_seas; NRMSE_ann_ann];
r_vals      = [R_ann_ra;       R_ann_seas;     R_ann_ann];
rd_vals     = [RD_ann_ra;      RD_ann_seas;    RD_ann_ann];

fig1d = figure('Position',[50 50 1200 480]);

subplot(1,3,1); hold on;
for k = 1:3, bar(k, nrmse_vals(k), 0.5, 'FaceColor', clrs(k,:)); end
set(gca,'XTick',1:3,'XTickLabel',methods,'FontSize',10);
ylabel('NRMSE','FontSize',11); title('NRMSE (annual)','FontWeight','bold');
grid on; box on;

subplot(1,3,2); hold on;
for k = 1:3, bar(k, r_vals(k), 0.5, 'FaceColor', clrs(k,:)); end
yline(1,'k--','LineWidth',1.2,'HandleVisibility','off');
set(gca,'XTick',1:3,'XTickLabel',methods,'FontSize',10);
ylim([0 1.1]); ylabel('Pearson r','FontSize',11); title('Pearson r (annual)','FontWeight','bold');
grid on; box on;

subplot(1,3,3); hold on;
for k = 1:3, bar(k, rd_vals(k), 0.5, 'FaceColor', clrs(k,:)); end
yline(0,'k-','LineWidth',1.5,'HandleVisibility','off');
set(gca,'XTick',1:3,'XTickLabel',methods,'FontSize',10);
ylabel('RD%','FontSize',11); title('RD% (annual)','FontWeight','bold');
grid on; box on;

sgtitle({sprintf('Annual Results (KED) -- CSV6  [SG1=%s, SG2=%s, Thr=%g]', ...
        mat2str(USE_SAFEGUARD_1), mat2str(USE_SAFEGUARD_2), SAFEGUARD_2_THRESHOLD), ...
    sprintf('NRMSE: Radar=%.3f | Seasonal KED=%.3f | Annual KED=%.3f', ...
        NRMSE_ann_ra, NRMSE_ann_seas, NRMSE_ann_ann)}, ...
    'FontWeight','bold','FontSize',11);

saveas(fig1d, fullfile(out_dir,'Fig_V1d_KED_Annual_Summary.png'));
savefig(fig1d, fullfile(out_dir,'Fig_V1d_KED_Annual_Summary.fig'));

%% ============================================================
%  نمودار ۱e: Bin-Averaged Scatter Plot
%% ============================================================
mask_sc = ~isnan(St_mat) & ~isnan(Ra_mat) & Ra_mat > 0 & ...
          ~isnan(Corr_seas) & ~isnan(Corr_ann);

Ra_sc   = Ra_mat(mask_sc);
Seas_sc = Corr_seas(mask_sc);
Ann_sc  = Corr_ann(mask_sc);
St_sc   = St_mat(mask_sc);

bin_edges = 0:0.2:60;
n_bins    = numel(bin_edges) - 1;

bin_x_ra   = nan(n_bins,1); bin_y_ra   = nan(n_bins,1); bin_n_ra   = zeros(n_bins,1);
bin_x_seas = nan(n_bins,1); bin_y_seas = nan(n_bins,1); bin_n_seas = zeros(n_bins,1);
bin_x_ann  = nan(n_bins,1); bin_y_ann  = nan(n_bins,1); bin_n_ann  = zeros(n_bins,1);

for b = 1:n_bins
    lo = bin_edges(b); hi = bin_edges(b+1);

    in_b = Ra_sc >= lo & Ra_sc < hi;
    if sum(in_b) >= 2
        bin_x_ra(b) = mean(Ra_sc(in_b)); bin_y_ra(b) = mean(St_sc(in_b)); bin_n_ra(b) = sum(in_b);
    end
    in_b = Seas_sc >= lo & Seas_sc < hi;
    if sum(in_b) >= 2
        bin_x_seas(b) = mean(Seas_sc(in_b)); bin_y_seas(b) = mean(St_sc(in_b)); bin_n_seas(b) = sum(in_b);
    end
    in_b = Ann_sc >= lo & Ann_sc < hi;
    if sum(in_b) >= 2
        bin_x_ann(b) = mean(Ann_sc(in_b)); bin_y_ann(b) = mean(St_sc(in_b)); bin_n_ann(b) = sum(in_b);
    end
end

ok_ra   = ~isnan(bin_x_ra)   & ~isnan(bin_y_ra);
ok_seas = ~isnan(bin_x_seas) & ~isnan(bin_y_seas);
ok_ann  = ~isnan(bin_x_ann)  & ~isnan(bin_y_ann);

tmp = corrcoef(bin_x_ra(ok_ra),     bin_y_ra(ok_ra));     r_bin_ra   = tmp(1,2);
tmp = corrcoef(bin_x_seas(ok_seas), bin_y_seas(ok_seas)); r_bin_seas = tmp(1,2);
tmp = corrcoef(bin_x_ann(ok_ann),   bin_y_ann(ok_ann));   r_bin_ann  = tmp(1,2);

sz_ra   = max(20, 15*log1p(bin_n_ra(ok_ra)));
sz_seas = max(20, 15*log1p(bin_n_seas(ok_seas)));
sz_ann  = max(20, 15*log1p(bin_n_ann(ok_ann)));

ax_max = max([ max(bin_x_ra(ok_ra),   [],'omitnan'); ...
               max(bin_x_seas(ok_seas),[],'omitnan'); ...
               max(bin_x_ann(ok_ann),  [],'omitnan'); ...
               max(bin_y_ra(ok_ra),    [],'omitnan'); ...
               max(bin_y_seas(ok_seas),[],'omitnan'); ...
               max(bin_y_ann(ok_ann),  [],'omitnan') ]);
ax_max = ceil(ax_max * 1.05);

fig1e = figure('Position',[50 50 750 650]); hold on;
h_ra   = scatter(bin_x_ra(ok_ra),     bin_y_ra(ok_ra),     sz_ra,   ...
    [0.6 0.6 0.6],    'filled', 'MarkerFaceAlpha',0.85, ...
    'DisplayName', sprintf('Before (r=%.3f)', r_bin_ra));
h_seas = scatter(bin_x_seas(ok_seas), bin_y_seas(ok_seas), sz_seas, ...
    [0.85 0.20 0.20], 'filled', 'MarkerFaceAlpha',0.85, ...
    'DisplayName', sprintf('Seasonal KED (r=%.3f)', r_bin_seas));
h_ann  = scatter(bin_x_ann(ok_ann),   bin_y_ann(ok_ann),   sz_ann,  ...
    [0.20 0.47 0.75], 'filled', 'MarkerFaceAlpha',0.85, ...
    'DisplayName', sprintf('Annual KED (r=%.3f)', r_bin_ann));

n_ra_ok   = bin_n_ra(ok_ra);
n_seas_ok = bin_n_seas(ok_seas);
n_ann_ok  = bin_n_ann(ok_ann);
rows_ra   = arrayfun(@(n) sprintf('N=%d hours', n), n_ra_ok,   'UniformOutput',false);
rows_seas = arrayfun(@(n) sprintf('N=%d hours', n), n_seas_ok, 'UniformOutput',false);
rows_ann  = arrayfun(@(n) sprintf('N=%d hours', n), n_ann_ok,  'UniformOutput',false);
h_ra.DataTipTemplate.DataTipRows(end+1)   = dataTipTextRow('', rows_ra);
h_seas.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('', rows_seas);
h_ann.DataTipTemplate.DataTipRows(end+1)  = dataTipTextRow('', rows_ann);

plot([0 ax_max],[0 ax_max],'k--','LineWidth',1.5,'DisplayName','1:1 line');

xlim([0 ax_max]); ylim([0 ax_max]);
xlabel('CPC bin mean (mm/h)','FontSize',11);
ylabel('Gauge mean (mm/h)','FontSize',11);
title({'[CSV6, KED] Bin-Avg Scatter -- Overlay [all GS, CPC>0]', ...
    sprintf('bin width=0.2 mm/h | Before r=%.3f | Seasonal KED r=%.3f | Annual KED r=%.3f', ...
        r_bin_ra, r_bin_seas, r_bin_ann)}, 'FontWeight','bold');
legend('Location','northwest','FontSize',9);
grid on; box on; axis square;

saveas(fig1e, fullfile(out_dir,'Fig_V1e_KED_BinAvg_Scatter.png'));
savefig(fig1e, fullfile(out_dir,'Fig_V1e_KED_BinAvg_Scatter.fig'));

%% ============================================================
%  نمودار ۲: کومولاتیو pooled
%% ============================================================
fig2 = figure('Position',[50 50 1000 480]); hold on;
plot(times_ok, cumSt,   '-',  'Color',[0 0 0],          'LineWidth',2.2, 'DisplayName','Station');
plot(times_ok, cumRa,   '--', 'Color',[0.6 0.6 0.6],    'LineWidth',1.5, 'DisplayName','Radar (before)');
plot(times_ok, cumSeas, '-',  'Color',[0.20 0.47 0.75], 'LineWidth',1.8, 'DisplayName','Seasonal KED (+fallback)');
plot(times_ok, cumAnn,  '--', 'Color',[0.47 0.67 0.19], 'LineWidth',1.8, 'DisplayName','Annual KED (+fallback)');
xlabel('Time'); ylabel('Cumulative spatial mean rainfall (mm)');
title('Cumulative hourly spatial mean (KED) -- CSV6 (extended grid, radar fallback)','FontWeight','bold');
legend('Location','northwest','FontSize',9); grid on; box on;
saveas(fig2, fullfile(out_dir,'Fig_V2_KED_CumulativeSum_Comparison.png'));
savefig(fig2, fullfile(out_dir,'Fig_V2_KED_CumulativeSum_Comparison.fig'));

%% ============================================================
%  نمودار ۳: تفاضل نسبی کومولاتیو pooled (از ۱۶ مارس)
%% ============================================================
start_plot = datetime(2025,3,16,0,0,0);
pm = times_ok >= start_plot;

fig3 = figure('Position',[50 50 1000 480]); hold on;
plot(times_ok(pm), relDiff_Ra(pm),   '--', 'Color',[0.6 0.6 0.6],    'LineWidth',1.5, ...
    'DisplayName',sprintf('Radar    (final: %.1f%%)', relDiff_Ra(end)));
plot(times_ok(pm), relDiff_Seas(pm), '-',  'Color',[0.20 0.47 0.75], 'LineWidth',2.0, ...
    'DisplayName',sprintf('Seasonal KED (final: %.1f%%)', relDiff_Seas(end)));
plot(times_ok(pm), relDiff_Ann(pm),  '--', 'Color',[0.47 0.67 0.19], 'LineWidth',2.0, ...
    'DisplayName',sprintf('Annual KED   (final: %.1f%%)', relDiff_Ann(end)));
yline(0,'k-','LineWidth',1.5,'HandleVisibility','off');
xlabel('Time'); ylabel('Cumulative relative difference (%)');
title({'Cumulative relative difference (KED) -- CSV6 (extended grid, radar fallback)', ...
    sprintf('from %s', datestr(start_plot,'dd mmm yyyy'))}, 'FontWeight','bold');
legend('Location','northeast','FontSize',9); grid on; box on;
saveas(fig3, fullfile(out_dir,'Fig_V3_KED_RelDiff_Comparison_from16Mar.png'));
savefig(fig3, fullfile(out_dir,'Fig_V3_KED_RelDiff_Comparison_from16Mar.fig'));

%% ============================================================
%  PER-STATION: کومولاتیو + تفاضل نسبی
%% ============================================================
cum_station = struct();

for v = 1:n6
    ok_v = ~isnan(St_mat(:,v)) & ~isnan(Ra_mat(:,v)) & ...
           ~isnan(Corr_seas(:,v)) & ~isnan(Corr_ann(:,v));
    times_v = grid_times(ok_v);

    cumSt_v   = cumsum(St_mat(ok_v,v));
    cumRa_v   = cumsum(Ra_mat(ok_v,v));
    cumSeas_v = cumsum(Corr_seas(ok_v,v));
    cumAnn_v  = cumsum(Corr_ann(ok_v,v));

    relDiff_Ra_v   = (cumRa_v   - cumSt_v) ./ cumSt_v * 100;
    relDiff_Seas_v = (cumSeas_v - cumSt_v) ./ cumSt_v * 100;
    relDiff_Ann_v  = (cumAnn_v  - cumSt_v) ./ cumSt_v * 100;

    cum_station(v).name         = csv6_names{v};
    cum_station(v).times        = times_v;
    cum_station(v).cumSt        = cumSt_v;
    cum_station(v).cumRa        = cumRa_v;
    cum_station(v).cumSeas      = cumSeas_v;
    cum_station(v).cumAnn       = cumAnn_v;
    cum_station(v).relDiff_Ra   = relDiff_Ra_v;
    cum_station(v).relDiff_Seas = relDiff_Seas_v;
    cum_station(v).relDiff_Ann  = relDiff_Ann_v;

    safe_name = matlab.lang.makeValidName(csv6_names{v});

    figV2_v = figure('Position',[50 50 1000 480]); hold on;
    plot(times_v, cumSt_v,   '-',  'Color',[0 0 0],          'LineWidth',2.2, 'DisplayName','Station');
    plot(times_v, cumRa_v,   '--', 'Color',[0.6 0.6 0.6],    'LineWidth',1.5, 'DisplayName','Radar (before)');
    plot(times_v, cumSeas_v, '-',  'Color',[0.20 0.47 0.75], 'LineWidth',1.8, 'DisplayName','Seasonal KED (+fallback)');
    plot(times_v, cumAnn_v,  '--', 'Color',[0.47 0.67 0.19], 'LineWidth',1.8, 'DisplayName','Annual KED (+fallback)');
    xlabel('Time'); ylabel('Cumulative rainfall (mm)');
    title(sprintf('Cumulative rainfall (KED) -- %s', csv6_names{v}), 'FontWeight','bold','Interpreter','none');
    legend('Location','northwest','FontSize',9); grid on; box on;
    saveas(figV2_v, fullfile(out_dir, sprintf('Fig_V2_KED_Station_%02d_%s.png', v, safe_name)));
    savefig(figV2_v, fullfile(out_dir, sprintf('Fig_V2_KED_Station_%02d_%s.fig', v, safe_name)));
    close(figV2_v);

    pm_v = times_v >= start_plot;
    figV3_v = figure('Position',[50 50 1000 480]); hold on;
    plot(times_v(pm_v), relDiff_Ra_v(pm_v),   '--', 'Color',[0.6 0.6 0.6],    'LineWidth',1.5, ...
        'DisplayName', sprintf('Radar    (final: %.1f%%)', relDiff_Ra_v(end)));
    plot(times_v(pm_v), relDiff_Seas_v(pm_v), '-',  'Color',[0.20 0.47 0.75], 'LineWidth',2.0, ...
        'DisplayName', sprintf('Seasonal KED (final: %.1f%%)', relDiff_Seas_v(end)));
    plot(times_v(pm_v), relDiff_Ann_v(pm_v),  '--', 'Color',[0.47 0.67 0.19], 'LineWidth',2.0, ...
        'DisplayName', sprintf('Annual KED   (final: %.1f%%)', relDiff_Ann_v(end)));
    yline(0,'k-','LineWidth',1.5,'HandleVisibility','off');
    xlabel('Time'); ylabel('Cumulative relative difference (%)');
    title(sprintf('Cumulative relative difference (KED) -- %s', csv6_names{v}), 'FontWeight','bold','Interpreter','none');
    legend('Location','northeast','FontSize',9); grid on; box on;
    saveas(figV3_v, fullfile(out_dir, sprintf('Fig_V3_KED_Station_%02d_%s.png', v, safe_name)));
    savefig(figV3_v, fullfile(out_dir, sprintf('Fig_V3_KED_Station_%02d_%s.fig', v, safe_name)));
    close(figV3_v);
end

%% ============================================================
%  PER-STATION NRMSE: ماهانه و سالانه
%% ============================================================
NRMSE_station_month_ra   = nan(12, n6);
NRMSE_station_month_seas = nan(12, n6);
NRMSE_station_month_ann  = nan(12, n6);
N_station_month          = zeros(12, n6);
NRMSE_station_annual_ra   = nan(n6,1);
NRMSE_station_annual_seas = nan(n6,1);
NRMSE_station_annual_ann  = nan(n6,1);
N_station_annual          = zeros(n6,1);

for v = 1:n6
    ok_v_all = ~isnan(St_mat(:,v)) & ~isnan(Ra_mat(:,v)) & Ra_mat(:,v) > 0 & ...
               ~isnan(Corr_seas(:,v)) & ~isnan(Corr_ann(:,v));

    N_station_annual(v) = sum(ok_v_all);
    if N_station_annual(v) >= 5
        St_v = St_mat(ok_v_all,v); Ra_v = Ra_mat(ok_v_all,v);
        Cs_v = Corr_seas(ok_v_all,v); Ca_v = Corr_ann(ok_v_all,v);
        nv_v = mean(St_v);
        NRMSE_station_annual_ra(v)   = sqrt(mean((Ra_v - St_v).^2)) / nv_v;
        NRMSE_station_annual_seas(v) = sqrt(mean((Cs_v - St_v).^2)) / nv_v;
        NRMSE_station_annual_ann(v)  = sqrt(mean((Ca_v - St_v).^2)) / nv_v;
    end
    for m = 1:12
        ok_vm = ok_v_all & (month_of_grid == m);
        N_station_month(m,v) = sum(ok_vm);
        if N_station_month(m,v) < 5, continue; end
        St_vm = St_mat(ok_vm,v); Ra_vm = Ra_mat(ok_vm,v);
        Cs_vm = Corr_seas(ok_vm,v); Ca_vm = Corr_ann(ok_vm,v);
        nv_vm = mean(St_vm);
        NRMSE_station_month_ra(m,v)   = sqrt(mean((Ra_vm - St_vm).^2)) / nv_vm;
        NRMSE_station_month_seas(m,v) = sqrt(mean((Cs_vm - St_vm).^2)) / nv_vm;
        NRMSE_station_month_ann(m,v)  = sqrt(mean((Ca_vm - St_vm).^2)) / nv_vm;
    end
end

%% ============================================================
%  ذخیره جداول
%% ============================================================
T_nrmse = table(month_names', N_month, NRMSE_ra, NRMSE_seas, NRMSE_ann, ...
    (NRMSE_ra-NRMSE_seas)./NRMSE_ra*100, (NRMSE_ra-NRMSE_ann)./NRMSE_ra*100, ...
    'VariableNames',{'Month','N_hours','NRMSE_Radar','NRMSE_Seasonal_KED','NRMSE_Annual_KED', ...
        'Improvement_Seasonal_pct','Improvement_Annual_pct'});
writetable(T_nrmse, fullfile(out_dir,'Table_NRMSE_Comparison_KED.csv'));

T_ann = table({'Radar';'Seasonal KED';'Annual KED'}, ...
    [NRMSE_ann_ra; NRMSE_ann_seas; NRMSE_ann_ann], ...
    [0; 100*(NRMSE_ann_ra-NRMSE_ann_seas)/NRMSE_ann_ra; 100*(NRMSE_ann_ra-NRMSE_ann_ann)/NRMSE_ann_ra], ...
    'VariableNames',{'Method','NRMSE_Annual','Improvement_pct'});
writetable(T_ann, fullfile(out_dir,'Table_NRMSE_Annual_Comparison_KED.csv'));

T_cum = table(times_ok, cumSt, cumRa, cumSeas, cumAnn, ...
    relDiff_Ra, relDiff_Seas, relDiff_Ann, ...
    'VariableNames',{'Time','CumSum_Station','CumSum_Radar', ...
        'CumSum_Seasonal_KED','CumSum_Annual_KED', ...
        'RelDiff_Radar_pct','RelDiff_Seasonal_KED_pct','RelDiff_Annual_KED_pct'});
writetable(T_cum, fullfile(out_dir,'Table_Cumulative_Comparison_KED.csv'));

station_long={}; time_long=[]; cumSt_long=[]; cumRa_long=[]; cumSeas_long=[]; cumAnn_long=[];
relRa_long=[]; relSeas_long=[]; relAnn_long=[];
for v = 1:n6
    nrows = numel(cum_station(v).times);
    station_long = [station_long; repmat({cum_station(v).name}, nrows, 1)]; %#ok<AGROW>
    time_long    = [time_long;    cum_station(v).times]; %#ok<AGROW>
    cumSt_long   = [cumSt_long;   cum_station(v).cumSt]; %#ok<AGROW>
    cumRa_long   = [cumRa_long;   cum_station(v).cumRa]; %#ok<AGROW>
    cumSeas_long = [cumSeas_long; cum_station(v).cumSeas]; %#ok<AGROW>
    cumAnn_long  = [cumAnn_long;  cum_station(v).cumAnn]; %#ok<AGROW>
    relRa_long   = [relRa_long;   cum_station(v).relDiff_Ra]; %#ok<AGROW>
    relSeas_long = [relSeas_long; cum_station(v).relDiff_Seas]; %#ok<AGROW>
    relAnn_long  = [relAnn_long;  cum_station(v).relDiff_Ann]; %#ok<AGROW>
end
T_cum_station = table(station_long, time_long, cumSt_long, cumRa_long, cumSeas_long, cumAnn_long, ...
    relRa_long, relSeas_long, relAnn_long, ...
    'VariableNames', {'Station','Time','CumSum_Station','CumSum_Radar', ...
        'CumSum_Seasonal_KED','CumSum_Annual_KED','RelDiff_Radar_pct','RelDiff_Seasonal_KED_pct','RelDiff_Annual_KED_pct'});
writetable(T_cum_station, fullfile(out_dir,'Table_Cumulative_PerStation_KED.csv'));

station_col={}; month_col={}; N_col=[]; ra_col=[]; seas_col=[]; ann_col=[];
for v = 1:n6
    for m = 1:12
        station_col{end+1,1} = csv6_names{v}; %#ok<AGROW>
        month_col{end+1,1}   = month_names{m}; %#ok<AGROW>
        N_col(end+1,1)       = N_station_month(m,v); %#ok<AGROW>
        ra_col(end+1,1)      = NRMSE_station_month_ra(m,v); %#ok<AGROW>
        seas_col(end+1,1)    = NRMSE_station_month_seas(m,v); %#ok<AGROW>
        ann_col(end+1,1)     = NRMSE_station_month_ann(m,v); %#ok<AGROW>
    end
end
T_station_month = table(station_col, month_col, N_col, ra_col, seas_col, ann_col, ...
    'VariableNames', {'Station','Month','N_hours','NRMSE_Radar','NRMSE_Seasonal_KED','NRMSE_Annual_KED'});
writetable(T_station_month, fullfile(out_dir,'Table_NRMSE_PerStation_PerMonth_KED.csv'));

T_station_annual = table(csv6_names(:), N_station_annual, ...
    NRMSE_station_annual_ra, NRMSE_station_annual_seas, NRMSE_station_annual_ann, ...
    100*(NRMSE_station_annual_ra-NRMSE_station_annual_seas)./NRMSE_station_annual_ra, ...
    100*(NRMSE_station_annual_ra-NRMSE_station_annual_ann)./NRMSE_station_annual_ra, ...
    'VariableNames', {'Station','N_hours','NRMSE_Radar','NRMSE_Seasonal_KED','NRMSE_Annual_KED', ...
                       'Improvement_Seasonal_pct','Improvement_Annual_pct'});
writetable(T_station_annual, fullfile(out_dir,'Table_NRMSE_PerStation_Annual_KED.csv'));

fprintf('\nAll KED outputs (CSVs + figures + tables) saved to: %s\n', out_dir);


%% ============================================================
%  تابع محلی: کریجینگ با درفت خارجی (KED) -- Safeguard قابل‌تنظیم
%  هدف = ln(S)  |  درفت = ln(R)  |  خروجی = exp(پیش‌بینی مستقیم S)
%% ============================================================

%% ============================================================
%  Local functions (moved to the end of the script so it also runs
%  on MATLAB releases older than R2024a)
%% ============================================================
function [corr_val, flag_extrap, max_lambda] = do_kriging_KED(t_idx, Sx_v, Sy_v, ...
        C0, CH, AH, rad, ...
        lnS_dense_tr, lnR_dense_tr, Dfull_tr, Sx_tr, Sy_tr, n_train, ...
        Ra_val, exp_model, min_stations, ...
        use_safeguard_1, use_safeguard_2, safeguard_2_threshold)

    corr_val = NaN; flag_extrap = false; max_lambda = NaN;

    d_to_train = sqrt((Sx_tr-Sx_v).^2 + (Sy_tr-Sy_v).^2);
    obs_mask = ~isnan(lnS_dense_tr(t_idx,:));
    if isinf(rad)
        in_rad = true(1,n_train);
    else
        in_rad = d_to_train <= rad;
    end
    use_idx = find(obs_mask & in_rad);
    if numel(use_idx) < min_stations, return; end

    no         = numel(use_idx);
    sill_total = C0 + CH;

    C_mat = zeros(no+2);
    for a = 1:no
        for b = 1:no
            h_ab = Dfull_tr(use_idx(a),use_idx(b));
            if h_ab == 0
                C_mat(a,b) = sill_total;
            else
                C_mat(a,b) = sill_total - exp_model(C0,CH,AH,h_ab);
            end
        end
    end

    drift_nb   = lnR_dense_tr(t_idx, use_idx)';
    target_lnR = log(Ra_val);
    flag_extrap = (target_lnR < min(drift_nb)) || (target_lnR > max(drift_nb));

    % --- Safeguard 1 (اختیاری): رد کردن برون‌یابی در درفت ---
    if use_safeguard_1 && flag_extrap
        return;
    end

    C_mat(no+1, 1:no) = 1;          C_mat(1:no, no+1) = 1;
    C_mat(no+2, 1:no) = drift_nb';  C_mat(1:no, no+2) = drift_nb;

    c_vec = zeros(no+2,1);
    for a = 1:no
        h_ia     = sqrt((Sx_tr(use_idx(a))-Sx_v)^2 + (Sy_tr(use_idx(a))-Sy_v)^2);
        c_vec(a) = sill_total - exp_model(C0,CH,AH,h_ia);
    end
    c_vec(no+1) = 1;
    c_vec(no+2) = target_lnR;

    if rcond(C_mat) < 1e-12, return; end
    lambda     = C_mat \ c_vec;
    max_lambda = max(abs(lambda));

    % --- Safeguard 2 (اختیاری): رد کردن وزن نجومی، آستانهٔ دلخواه ---
    if use_safeguard_2 && max_lambda > safeguard_2_threshold
        return;
    end

    lnS_pred = lambda(1:no)' * lnS_dense_tr(t_idx,use_idx)';
    corr_val = exp(lnS_pred);
end
