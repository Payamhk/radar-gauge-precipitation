%% ------------------------------------------------------------------
%  qc_hourly_consistency_check.m
%  Part of: radar-gauge-precipitation (Payam Heidarian, University of
%  Brescia / MeteoSwiss collaboration). Inline comments are partly in
%  Persian. Data are not included; set the <DATA_ROOT> paths below to
%  your own folders before running.
%% ------------------------------------------------------------------
%% ============================================================
%  Consorzio Hourly QC - v2 FINAL + Q95 VARIANT (4 کلاس کاملاً جدا: NoRain / Light / Moderate / Heavy)
%
%  این نسخه دقیقاً همان منطق فایل اصلی (v2 FINAL) را دارد، با یک اضافه:
%  علاوه بر Consorzio_hourly_Final_Table (که در آن Light/Moderate/Heavy
%  اصلاً بر اساس SARD فلگ نمی‌شوند)، یک جدول دوم به نام
%  Consorzio_hourly_Final_Table_95 هم ساخته می‌شود که در آن رکوردهای
%  Light/Moderate/Heavy با SARD (metric قدرمطلق نسبی) >= صدک ۹۵ درصد
%  (محاسبه‌شده مستقل برای هر کلاس، از روی pool همان کلاس) به عنوان
%  NotAligned فلگ می‌شوند و از Valid3 آن جدول حذف می‌شوند.
%
%  قانون NoRain (آستانه‌ی ثابت ۳mm روی P_IDW) و قانون تکمیلی
%  AllNeighboursZero (آستانه‌ی ثابت ۱mm) در هر دو جدول یکسان و بدون
%  تغییر باقی می‌مانند.
%
%  Progressive filtering (ترتیب توافق‌شده):
%
%    Total (raw aggregated hourly rain, artificial/missed hours kept as 0
%           and treated as ordinary valid values everywhere)
%      |
%      v
%    Step 1: FlagRepeat4h  = 4+ consecutive equal NONZERO hours, checked on Total
%    Valid1 = Total - FlagRepeat4h
%      |
%      v
%    Step 2: Flag200 = P > Hourly_threshold (200mm/h), checked on Valid1
%    Valid2 = Valid1 - Flag200
%      |
%      v
%    Step 3: Consistency check, checked on Valid2. Neighbours are also
%             filtered to Valid2 before being used in IDW.
%
%       - If number of VALID2 neighbours < 3  -> NotEvaluated (excluded from
%         Valid3, but NOT flagged as an error; a third, neutral status)
%       - If number of VALID2 neighbours >= 3 -> classify target rain (WMO)
%             into FOUR SEPARATE classes:
%
%             NoRain    (P = 0)          -> metric = P_IDW
%                 FIXED threshold: flag if P_IDW >= NoRain_fixed_threshold_mm (3 mm)
%             Light     (0 < P < 2.5)    -> metric = |2*(P_target-P_IDW)/(P_target+P_IDW)|  (SARD)
%             Moderate  (2.5 <= P < 10)  -> metric = |2*(P_target-P_IDW)/(P_target+P_IDW)|  (SARD)
%             Heavy     (P >= 10)        -> metric = |2*(P_target-P_IDW)/(P_target+P_IDW)|  (SARD)
%
%    Valid3 (final, ORIGINAL table)    = Valid2 - NotEvaluated - FlagNotAligned            (L/M/H NOT flagged by SARD)
%    Valid3_95 (final, Q95 table)      = Valid2 - NotEvaluated - FlagNotAligned_95          (L/M/H flagged if SARD >= class-specific Q95)
%
%  Artificial (missed) hours: informational column only. Their value (0) is
%  never removed from Valid1/Valid2/Valid3 and is never treated as an outlier.
%% ============================================================
clear, clc

tic

%% ---- 0. Parameters ----
start_datetime      = datetime(2025,03,13,0,0,0);
end_time             = datetime(2026,03,13,0,0,0);
Hourly_threshold     = 200;     % mm/h  (Step 2)
Repeat_run_length    = 4;       % consecutive equal nonzero hours (Step 1)
Min_valid_neighbours = 3;       % below this -> NotEvaluated (Step 3)

% WMO intensity class breakpoints (4 کلاس مستقل)
Light_upper_mm    = 2.5;    % 0 < P < 2.5      -> Light
Moderate_upper_mm = 10.0;   % 2.5 <= P < 10     -> Moderate
% P >= 10 -> Heavy، P == 0 -> No Rain

% ---- No Rain class: آستانه‌ی ثابت (نه صدک‌محور) ----
NoRain_fixed_threshold_mm = 3.0;

% ---- قانون مکمل (جهت معکوس No-Rain): همه‌ی همسایه‌های معتبر صفر ----
AllNeighboursZero_target_threshold_mm = 1.0;

% ---- NEW: صدک اعمال‌شده برای Light/Moderate/Heavy در جدول دوم ----
Q95_level = 0.95;

Consersio_coordinates = readtable("<DATA_ROOT>\Moghayese gs vs nc franciacorta\چیزی که استاد خواسته\ایستگاه های زمینی دادم به استاد\16-3-2026 یکساله\44 Stations برای یکسال.csv");

Consorzio_Raw_Controlled = struct2table(load('Output_Consorzio_Raw/Consorzio_Final_Output_Raw44sbyProf.mat'));
Consorzio_Raw_Controlled = Consorzio_Raw_Controlled{:,1};

n_stations = height(Consorzio_Raw_Controlled{:,1});

%% ---- 1. ساخت گرید زمانی ساعتی ----
ORDINAL_TIMES_Hourly = (start_datetime : hours(1) : end_time)';
time_edges           = ORDINAL_TIMES_Hourly + seconds(1);
time_centers         = ORDINAL_TIMES_Hourly(2:end);
n_hours              = numel(time_centers);

%% ---- 2. Aggregate به ساعت (Total) + شناسایی Artificial ----
Consorzio_hourly_Final = cell(n_stations, 2);
artificial_matrix      = false(n_hours, n_stations);
rain_matrix             = zeros(n_hours, n_stations);

for i = 1:n_stations
    fprintf('Hourly aggregation: station %d/%d\n', i, n_stations);

    Consorzio_station_name = Consorzio_Raw_Controlled{i,1};
    raw_sub  = Consorzio_Raw_Controlled{i,2}{1};
    raw_time = raw_sub{:,1};
    raw_rain = raw_sub{:,2};
    raw_qual = raw_sub{:,3};

    [~, ~, bin] = histcounts(datenum(raw_time), datenum(time_edges));
    valid_bin = bin >= 1 & bin <= n_hours;
    bin       = bin(valid_bin);
    raw_rain  = raw_rain(valid_bin);
    raw_qual  = raw_qual(valid_bin);

    rain_hourly = accumarray(bin, raw_rain, [n_hours 1], @sum, 0);

    is_missed        = contains(raw_qual, "missed time and precipitation");
    n_total_per_bin  = accumarray(bin, ones(size(bin)), [n_hours 1], @sum, 0);
    n_missed_per_bin = accumarray(bin, double(is_missed), [n_hours 1], @sum, 0);
    all_missed       = (n_total_per_bin > 0) & (n_missed_per_bin == n_total_per_bin);

    rain_hourly(all_missed) = 0;

    Consorzio_hourly        = table();
    Consorzio_hourly.time   = time_centers;
    Consorzio_hourly.rain   = rain_hourly;

    Consorzio_hourly_Final{i,1} = Consorzio_station_name;
    Consorzio_hourly_Final{i,2} = Consorzio_hourly;

    rain_matrix(:,i)        = rain_hourly;
    artificial_matrix(:,i)  = all_missed;
end

%% ---- 3. پیدا کردن 4 همسایه ثابت (نزدیک‌ترین بر اساس فاصله مکانی) ----
fprintf('Computing nearest neighbours...\n');
Consersio_coordinates.NearestStations = cell(height(Consersio_coordinates), 1);

for i = 1:height(Consersio_coordinates)
    d    = sqrt((Consersio_coordinates{:,2} - Consersio_coordinates{i,2}).^2 + ...
                (Consersio_coordinates{:,3} - Consersio_coordinates{i,3}).^2);
    d(i) = inf;
    [ds, idx] = mink(d, 5);
    Consersio_coordinates.NearestStations{i} = table( ...
        Consersio_coordinates{idx,1}, ds(:), ...
        Consersio_coordinates{idx,2}, ...
        Consersio_coordinates{idx,3}, ...
        'VariableNames', {'Station','Distance','X','Y'});
end

warning off
nb_col_all  = cell(n_stations, 1);
nb_dist_all = cell(n_stations, 1);
nb_xy_all   = cell(n_stations, 1);   % X,Y مستقیم از همون جدول NearestStations - بدون تطبیق دوباره

for i = 1:n_stations
    nb_table = Consersio_coordinates.NearestStations{i};
    nb_names = nb_table{:,1};
    nb_dist  = nb_table{:,2};
    nb_xy    = nb_table{:,3:4};   % X, Y مستقیماً از همین جدول
    n_nb     = numel(nb_names);

    nb_col = zeros(n_nb, 1);
    for n = 1:n_nb
        for ii = 1:n_stations
            if ii == i, continue; end
            if contains(lower(regexprep(string(Consorzio_hourly_Final{ii,1}),'[ _]','')), ...
                        lower(regexprep(string(nb_names(n)),'[ _]','')))
                nb_col(n) = ii;
                break;
            end
        end
    end

    found   = nb_col > 0;
    nb_col  = nb_col(found);
    nb_dist = nb_dist(found);
    nb_xy   = nb_xy(found, :);
    nb_col  = nb_col(1:min(4, end));
    nb_dist = nb_dist(1:min(4, end));
    nb_xy   = nb_xy(1:min(4, end), :);

    nb_col_all{i}  = nb_col;
    nb_dist_all{i} = nb_dist;
    nb_xy_all{i}   = nb_xy;
end

%% ---- 4. Step 1: FlagRepeat4h روی Total ----
fprintf('Step 1: repeated precipitation check (on Total)...\n');

FlagRepeat4h_matrix = false(n_hours, n_stations);

for i = 1:n_stations
    rain_series = rain_matrix(:,i);

    is_same    = [false; diff(rain_series) == 0] & (rain_series ~= 0);
    run_id     = cumsum(~is_same);
    run_len    = accumarray(run_id, ones(n_hours,1));
    len_series = run_len(run_id);

    FlagRepeat4h_matrix(:,i) = (len_series >= Repeat_run_length) & (rain_series ~= 0);
end

Valid1_matrix = ~FlagRepeat4h_matrix;   % Valid1 = Total - FlagRepeat4h

%% ---- 5. Step 2: Flag200 روی Valid1 ----
fprintf('Step 2: >%dmm/h check (on Valid1)...\n', Hourly_threshold);

Flag200_matrix = Valid1_matrix & (rain_matrix > Hourly_threshold);
Valid2_matrix  = Valid1_matrix & ~Flag200_matrix;   % Valid2 = Valid1 - Flag200

%% ---- 6. Step 3: Consistency check روی Valid2 (4 کلاس کاملاً مستقل) ----
fprintf('Step 3: consistency check (IDW / 4 WMO classes) on Valid2...\n');

NumValidNb_mat       = zeros(n_hours, n_stations);
PIDW_mat             = nan(n_hours, n_stations);
Eligible_matrix      = false(n_hours, n_stations);
NotEvaluated_matrix  = false(n_hours, n_stations);

Metric_NoRain_mat   = nan(n_hours, n_stations);
Metric_Light_mat    = nan(n_hours, n_stations);
Metric_Moderate_mat = nan(n_hours, n_stations);
Metric_Heavy_mat    = nan(n_hours, n_stations);
Class_mat           = strings(n_hours, n_stations);

% قانون مکمل: همه‌ی همسایه‌های معتبر صفر ولی هدف >= آستانه (مستقل از کلاس/SARD)
FlagAllNeighboursZero_matrix = false(n_hours, n_stations);

for i = 1:n_stations
    nb_col  = nb_col_all{i};
    nb_dist = nb_dist_all{i};
    n_nb    = numel(nb_col);

    if n_nb == 0
        NotEvaluated_matrix(:,i) = Valid2_matrix(:,i);
        continue;
    end

    nb_rain   = rain_matrix(:, nb_col);
    nb_valid2 = Valid2_matrix(:, nb_col);
    num_valid = sum(nb_valid2, 2);

    w       = 1 ./ (nb_dist .^ 2);
    w_mat   = repmat(w', n_hours, 1);
    w_valid = w_mat .* nb_valid2;
    sum_w   = sum(w_valid, 2);
    P_IDW   = sum(w_valid .* nb_rain, 2) ./ max(sum_w, eps);

    rain_i = rain_matrix(:,i);

    is_v2     = Valid2_matrix(:,i);
    enough_nb = num_valid >= Min_valid_neighbours;

    eligible_i = is_v2 & enough_nb;
    noteval_i  = is_v2 & ~enough_nb;

    NumValidNb_mat(:,i)      = num_valid;
    PIDW_mat(:,i)            = P_IDW;
    Eligible_matrix(:,i)     = eligible_i;
    NotEvaluated_matrix(:,i) = noteval_i;

    % ---- 4 کلاس کاملاً مستقل ----
    class_norain   = eligible_i & (rain_i == 0);
    class_light    = eligible_i & (rain_i > 0 & rain_i < Light_upper_mm);
    class_moderate = eligible_i & (rain_i >= Light_upper_mm & rain_i < Moderate_upper_mm);
    class_heavy    = eligible_i & (rain_i >= Moderate_upper_mm);

    Class_mat(class_norain,   i) = "NoRain";
    Class_mat(class_light,    i) = "Light";
    Class_mat(class_moderate, i) = "Moderate";
    Class_mat(class_heavy,    i) = "Heavy";

    Metric_NoRain_mat(class_norain, i)  = P_IDW(class_norain);

    % ---- Light/Moderate/Heavy: هر سه از یک فرمول واحد (SARD = symmetric absolute relative difference) ----
    % نکته: وقتی P_IDW=0 باشد (ولی P_target>0)، RD همیشه دقیقاً ریاضی 2
    % می‌شود (سقف فرمول)، صرف‌نظر از مقدار واقعی P_target. این رکوردها
    % ("isolated non-zero") از محاسبه‌ی Metric حذف می‌شوند تا هیستوگرام/
    % صدک را با یک مقدار ثابت و بی‌معنی خراب نکنند؛ به‌جای آن، این حالت
    % توسط قانون AllNeighboursZero (پایین‌تر) با آستانه‌ی مستقل خودش
    % بررسی می‌شود.
    class_light_valid    = class_light    & (P_IDW ~= 0);
    class_moderate_valid = class_moderate & (P_IDW ~= 0);
    class_heavy_valid    = class_heavy    & (P_IDW ~= 0);

    Metric_Light_mat(class_light_valid, i)    = abs( 2 * (rain_i(class_light_valid) - P_IDW(class_light_valid)) ...
                                            ./ max(rain_i(class_light_valid) + P_IDW(class_light_valid), eps) );
    Metric_Moderate_mat(class_moderate_valid, i) = abs( 2 * (rain_i(class_moderate_valid) - P_IDW(class_moderate_valid)) ...
                                            ./ max(rain_i(class_moderate_valid) + P_IDW(class_moderate_valid), eps) );
    Metric_Heavy_mat(class_heavy_valid, i) = abs( 2 * (rain_i(class_heavy_valid) - P_IDW(class_heavy_valid)) ...
                                            ./ max(rain_i(class_heavy_valid) + P_IDW(class_heavy_valid), eps) );

    % ---- قانون مکمل: همه‌ی همسایه‌های معتبر دقیقاً صفر، ولی هدف >= آستانه ----
    nb_nonzero_count = sum(nb_valid2 & (nb_rain > 0), 2);
    FlagAllNeighboursZero_matrix(:,i) = eligible_i & ...
        (rain_i >= AllNeighboursZero_target_threshold_mm) & ...
        (nb_nonzero_count == 0);
end

% ---- Pool های مستقل هر 4 کلاس ----
pool_norain   = Metric_NoRain_mat(~isnan(Metric_NoRain_mat));      % فقط برای گزارش/نمودار
pool_light    = Metric_Light_mat(~isnan(Metric_Light_mat));
pool_moderate = Metric_Moderate_mat(~isnan(Metric_Moderate_mat));
pool_heavy    = Metric_Heavy_mat(~isnan(Metric_Heavy_mat));

% ---- No Rain: آستانه ثابت (نه صدک‌محور) ----
thr_norain_upper = NoRain_fixed_threshold_mm;

fprintf('  No Rain  : n=%d, FIXED threshold (P_IDW) = %.3f mm  [not data-driven]\n', ...
    numel(pool_norain), thr_norain_upper);
fprintf('  Light    : n=%d  (no single threshold in ORIGINAL table; Q95 used in _95 table)\n', numel(pool_light));
fprintf('  Moderate : n=%d  (no single threshold in ORIGINAL table; Q95 used in _95 table)\n', numel(pool_moderate));
fprintf('  Heavy    : n=%d  (no single threshold in ORIGINAL table; Q95 used in _95 table)\n', numel(pool_heavy));

%% ---- 6a. آستانه‌های Q95 برای Light/Moderate/Heavy (فقط برای جدول _95) ----
thr_light_q95    = quantile(pool_light,    Q95_level);
thr_moderate_q95 = quantile(pool_moderate, Q95_level);
thr_heavy_q95    = quantile(pool_heavy,    Q95_level);

fprintf('\n  ---- Q95 thresholds on SARD (NEW, used ONLY in Consorzio_hourly_Final_Table_95) ----\n');
fprintf('  Light    Q95 (SARD) = %.4f\n', thr_light_q95);
fprintf('  Moderate Q95 (SARD) = %.4f\n', thr_moderate_q95);
fprintf('  Heavy    Q95 (SARD) = %.4f\n', thr_heavy_q95);

%% ---- 6b. ORIGINAL FlagNotAligned_matrix (Light/Moderate/Heavy NOT flagged by SARD) ----
FlagNotAligned_matrix = false(n_hours, n_stations);

mask_norain   = (Class_mat == "NoRain");
mask_light    = (Class_mat == "Light");
mask_moderate = (Class_mat == "Moderate");
mask_heavy    = (Class_mat == "Heavy");

FlagNotAligned_matrix(mask_norain) = Metric_NoRain_mat(mask_norain) >= thr_norain_upper;

FlagNotAligned_matrix = FlagNotAligned_matrix & Eligible_matrix;

% ---- اعمال قانون مکمل (OR با نتیجه‌ی قبلی، نه جایگزین) ----
FlagNotAligned_matrix = FlagNotAligned_matrix | FlagAllNeighboursZero_matrix;

fprintf('\n  [ORIGINAL] All-neighbours-zero rule: %d records flagged (target >= %.1fmm, all valid neighbours = 0)\n', ...
    sum(FlagAllNeighboursZero_matrix(:)), AllNeighboursZero_target_threshold_mm);
fprintf('  [ORIGINAL] Total flagged (NotAligned): %d records\n', sum(FlagNotAligned_matrix(:)));

%% ---- 6c. NEW: FlagNotAligned_matrix_95 (Light/Moderate/Heavy flagged if SARD >= class Q95) ----
FlagNotAligned_matrix_95 = false(n_hours, n_stations);

% NoRain: همان قانون ثابت قبلی (بدون تغییر)
FlagNotAligned_matrix_95(mask_norain) = Metric_NoRain_mat(mask_norain) >= thr_norain_upper;

% Light/Moderate/Heavy: این‌بار با آستانه‌ی Q95 روی SARD (مستقل برای هر کلاس)
FlagNotAligned_matrix_95(mask_light)    = Metric_Light_mat(mask_light)       >= thr_light_q95;
FlagNotAligned_matrix_95(mask_moderate) = Metric_Moderate_mat(mask_moderate) >= thr_moderate_q95;
FlagNotAligned_matrix_95(mask_heavy)    = Metric_Heavy_mat(mask_heavy)       >= thr_heavy_q95;

FlagNotAligned_matrix_95 = FlagNotAligned_matrix_95 & Eligible_matrix;

% قانون مکمل (همان قبلی، مستقل از کلاس) هم در نسخه‌ی Q95 اعمال می‌شود
FlagNotAligned_matrix_95 = FlagNotAligned_matrix_95 | FlagAllNeighboursZero_matrix;

fprintf('\n  [Q95] Total flagged (NotAligned_95): %d records\n', sum(FlagNotAligned_matrix_95(:)));
fprintf('    of which NoRain rule       : %d\n', sum(mask_norain & FlagNotAligned_matrix_95));
fprintf('    of which Light SARD>=Q95   : %d\n', sum(mask_light    & (Metric_Light_mat    >= thr_light_q95)    & Eligible_matrix));
fprintf('    of which Moderate SARD>=Q95: %d\n', sum(mask_moderate & (Metric_Moderate_mat >= thr_moderate_q95) & Eligible_matrix));
fprintf('    of which Heavy SARD>=Q95   : %d\n', sum(mask_heavy    & (Metric_Heavy_mat    >= thr_heavy_q95)    & Eligible_matrix));
fprintf('    of which AllNeighboursZero : %d\n', sum(FlagAllNeighboursZero_matrix(:)));

%% ---- 7. Valid3 نهایی (هر دو نسخه) ----
Valid3_matrix    = Valid2_matrix & ~NotEvaluated_matrix & ~FlagNotAligned_matrix;      % ORIGINAL
Valid3_matrix_95 = Valid2_matrix & ~NotEvaluated_matrix & ~FlagNotAligned_matrix_95;   % Q95 variant

%% ---- 7a. ماتریس Metric ترکیبی (هر رکورد مقدار Metric مخصوص کلاس خودش را می‌گیرد) ----
% از این هم در پلات‌ها و هم به‌عنوان ستون مرجع در جداول خروجی استفاده می‌شود.
Metric_Combined_mat = nan(n_hours, n_stations);
Metric_Combined_mat(Class_mat == "NoRain")   = Metric_NoRain_mat(Class_mat == "NoRain");
Metric_Combined_mat(Class_mat == "Light")    = Metric_Light_mat(Class_mat == "Light");
Metric_Combined_mat(Class_mat == "Moderate") = Metric_Moderate_mat(Class_mat == "Moderate");
Metric_Combined_mat(Class_mat == "Heavy")    = Metric_Heavy_mat(Class_mat == "Heavy");

%% ---- 7b. Diagnostic plots (4 کلاس کاملاً جدا، بر اساس نسخه‌ی ORIGINAL) ----
fprintf('Generating QC diagnostic plots (Step 3, 4 classes, Q90/Q95/Q99 scatters)...\n');

plot_folder = fullfile('Output_Consorzio_hourly_v3', 'QC_Diagnostic_Plots');
if ~exist(plot_folder, 'dir'), mkdir(plot_folder); end

all_class   = Class_mat(:);
all_rain    = rain_matrix(:);
all_pidw    = PIDW_mat(:);
all_flagged = FlagNotAligned_matrix(:);     % نمودارهای تشخیصی بر اساس نسخه‌ی ORIGINAL رنگ‌بندی می‌شوند
all_metric  = Metric_Combined_mat(:);

% سطوح کوانتایل مورد نظر برای اسکترهای جدا
quantile_levels = [0.90, 0.95, 0.99];

% چهار کلاس کاملاً مستقل - بدون هیچ ادغامی
classes_info = struct( ...
    'name',         {"NoRain",             "Light",                      "Moderate",                   "Heavy"}, ...
    'label',        {"No Rain (P = 0)",     "Light (0 < P < 2.5 mm)",     "Moderate (2.5 <= P < 10 mm)","Heavy (P >= 10 mm)"}, ...
    'metric_pool',  {pool_norain,           pool_light,                   pool_moderate,                pool_heavy}, ...
    'thr_low',      {NaN,                   NaN,                          NaN,                          NaN}, ...
    'thr_up',       {thr_norain_upper,      NaN,                          NaN,                          NaN}, ...
    'metric_label', {"P_{IDW} (mm)",        "SARD", "SARD", "SARD"}, ...
    'one_sided',    {true,                  true,                         true,                         true});

quantile_colors = [0.95 0.65 0.00; 0.90 0.10 0.10; 0.30 0.00 0.55];   % طلایی/کهربایی - قرمز روشن - بنفش تیره (کنتراست بیشتر بین سه رنگ)

% ---- تنظیمات فونت مشترک برای همه‌ی نمودارهای تشخیصی (بدون بولد، ۴ شماره بزرگ‌تر) ----
axis_fontsize   = 14;   % اعداد و تیک‌های محورها
label_fontsize  = 14;   % xlabel / ylabel
title_fontsize  = 15;   % عنوان نمودار
legend_fontsize = 13;   % لجند
annot_fontsize  = 13;   % متن‌های داخل نمودار (خطوط کوانتایل، QC threshold)

for c = 1:numel(classes_info)
    cname        = classes_info(c).name;
    clabel       = classes_info(c).label;
    metric_pool  = classes_info(c).metric_pool;
    thr_low      = classes_info(c).thr_low;
    thr_up       = classes_info(c).thr_up;
    metric_label = classes_info(c).metric_label;
    one_sided    = classes_info(c).one_sided;

    if isempty(metric_pool), continue; end

    idx_class = find(all_class == cname & ~isnan(all_metric));
    if isempty(idx_class), continue; end

    x_pidw   = all_pidw(idx_class);
    y_rain   = all_rain(idx_class);
    is_flag  = all_flagged(idx_class);
    m_class  = all_metric(idx_class);

    % ---- محاسبه‌ی آستانه‌های Q90/Q95/Q99 (برای هر کلاس، مستقل از threshold واقعی QC) ----
    q_thr_up  = nan(1, numel(quantile_levels));
    q_thr_low = nan(1, numel(quantile_levels));
    for qi = 1:numel(quantile_levels)
        p = quantile_levels(qi);
        q_thr_up(qi) = quantile(metric_pool, p);
        if ~one_sided
            q_thr_low(qi) = quantile(metric_pool, 1-p);
        end
    end

    % ================= هیستوگرام + خطوط کوانتایل =================
    fig1 = figure('Visible','off','Position',[100 100 900 550]);
    histogram(metric_pool, 60, 'FaceColor',[0.20 0.45 0.70], 'EdgeColor','none');
    hold on
    for qi = 1:numel(quantile_levels)
        clr = quantile_colors(qi,:);
        xline(q_thr_up(qi), '--', 'Color', clr, 'LineWidth',1.4);
        if ~one_sided
            xline(q_thr_low(qi), '--', 'Color', clr, 'LineWidth',1.4);
        end
    end
    yl = ylim;
    yMid = yl(1) + 0.5*(yl(2)-yl(1));
    yStep = 0.05*(yl(2)-yl(1));
    nQ = numel(quantile_levels);
    for qi = 1:nQ
        clr = quantile_colors(qi,:);
        yPos = yMid + yStep*(qi - (nQ+1)/2);
        text(q_thr_up(qi), yPos, sprintf('  Q%g = %.3f', quantile_levels(qi)*100, q_thr_up(qi)), ...
            'Color', clr, 'FontWeight','normal', 'FontSize', annot_fontsize, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
        if ~one_sided
            text(q_thr_low(qi), yPos, sprintf('Q%g = %.3f  ', (1-quantile_levels(qi))*100, q_thr_low(qi)), ...
                'Color', clr, 'FontWeight','normal', 'FontSize', annot_fontsize, 'HorizontalAlignment','right', 'VerticalAlignment','middle');
        end
    end
    if one_sided && ~isnan(thr_up)
        switch cname
            case "NoRain"
                qc_label = sprintf('QC thr (FIXED) = %.3f mm', thr_up);
            otherwise
                qc_label = sprintf('QC thr (|RD|) = %.3f', thr_up);
        end
        xline(thr_up, 'k-', qc_label, ...
            'LineWidth',2, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'FontWeight','normal', 'FontSize', annot_fontsize);
    elseif ~one_sided
        xline(thr_low, 'k-', sprintf('QC thr lower = %.3f', thr_low), ...
            'LineWidth',2, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'FontWeight','normal', 'FontSize', annot_fontsize);
        xline(thr_up, 'k-', sprintf('QC thr upper = %.3f', thr_up), ...
            'LineWidth',2, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'FontWeight','normal', 'FontSize', annot_fontsize);
    end
    title(sprintf('Consistency metric distribution - %s (n = %d)', clabel, numel(metric_pool)), 'FontWeight','normal', 'FontSize', title_fontsize);
    xlabel(metric_label, 'FontWeight','normal', 'FontSize', label_fontsize); ylabel('Frequency', 'FontWeight','normal', 'FontSize', label_fontsize);
    set(gca, 'FontWeight','normal', 'FontSize', axis_fontsize);
    grid on; box on
    saveas(fig1, fullfile(plot_folder, sprintf('Histogram_%s.png', cname)));
    fig1.Visible = 'on';
    savefig(fig1, fullfile(plot_folder, sprintf('Histogram_%s.fig', cname)));
    close(fig1);

    % ================= CDF =================
    fig1b = figure('Visible','off','Position',[100 100 900 550]);
    [f_cdf, x_cdf] = ecdf(metric_pool);
    plot(x_cdf, f_cdf*100, 'LineWidth',1.8, 'Color',[0.20 0.45 0.70]);
    hold on
    for qi = 1:numel(quantile_levels)
        clr = quantile_colors(qi,:);
        xline(q_thr_up(qi), '--', 'Color', clr, 'LineWidth',1.2);
        if ~one_sided
            xline(q_thr_low(qi), '--', 'Color', clr, 'LineWidth',1.2);
        end
    end
    yl = ylim;
    yMid = yl(1) + 0.5*(yl(2)-yl(1));
    yStep = 0.05*(yl(2)-yl(1));
    nQ = numel(quantile_levels);
    for qi = 1:nQ
        clr = quantile_colors(qi,:);
        yPos = yMid + yStep*(qi - (nQ+1)/2);
        text(q_thr_up(qi), yPos, sprintf('  Q%g = %.3f', quantile_levels(qi)*100, q_thr_up(qi)), ...
            'Color', clr, 'FontWeight','normal', 'FontSize', annot_fontsize, 'HorizontalAlignment','left', 'VerticalAlignment','middle');
        if ~one_sided
            text(q_thr_low(qi), yPos, sprintf('Q%g = %.3f  ', (1-quantile_levels(qi))*100, q_thr_low(qi)), ...
                'Color', clr, 'FontWeight','normal', 'FontSize', annot_fontsize, 'HorizontalAlignment','right', 'VerticalAlignment','middle');
        end
    end
    if one_sided && ~isnan(thr_up)
        switch cname
            case "NoRain"
                qc_label_cdf = sprintf('QC thr (FIXED) = %.3f mm', thr_up);
            otherwise
                qc_label_cdf = sprintf('QC thr (|RD|) = %.3f', thr_up);
        end
        xline(thr_up, 'k-', qc_label_cdf, ...
            'LineWidth',1.8, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'FontWeight','normal', 'FontSize', annot_fontsize);
    elseif ~one_sided
        xline(thr_low, 'k-', sprintf('QC thr lower = %.3f', thr_low), ...
            'LineWidth',1.8, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'FontWeight','normal', 'FontSize', annot_fontsize);
        xline(thr_up, 'k-', sprintf('QC thr upper = %.3f', thr_up), ...
            'LineWidth',1.8, 'LabelVerticalAlignment','middle', 'LabelOrientation','horizontal', 'FontWeight','normal', 'FontSize', annot_fontsize);
    end
    title(sprintf('Consistency metric CDF - %s (n = %d)', clabel, numel(metric_pool)), 'FontWeight','normal', 'FontSize', title_fontsize);
    xlabel(metric_label, 'FontWeight','normal', 'FontSize', label_fontsize); ylabel('Cumulative %', 'FontWeight','normal', 'FontSize', label_fontsize);
    set(gca, 'FontWeight','normal', 'FontSize', axis_fontsize);
    grid on; box on
    saveas(fig1b, fullfile(plot_folder, sprintf('CDF_%s.png', cname)));
    fig1b.Visible = 'on';
    savefig(fig1b, fullfile(plot_folder, sprintf('CDF_%s.fig', cname)));
    close(fig1b);

    % ================= اسکتر کلی (همه‌ی رکوردها، رنگ بر اساس Flag ORIGINAL) =================
    fig2 = figure('Visible','off','Position',[100 100 800 750]);
    hold on
    scatter(x_pidw(~is_flag), y_rain(~is_flag), 28, [0.20 0.50 0.80], 'filled', ...
        'MarkerFaceAlpha',0.35, 'DisplayName','Consistent (in Valid3)');
    scatter(x_pidw(is_flag),  y_rain(is_flag),  50, [0.85 0.10 0.10], 'filled', ...
        'MarkerFaceAlpha',0.85, 'DisplayName','Flagged (NotAligned)');
    max_ax = max([x_pidw; y_rain], [], 'omitnan') * 1.05;
    if isempty(max_ax) || max_ax <= 0 || isnan(max_ax), max_ax = 1; end
    plot([0 max_ax], [0 max_ax], 'k--', 'LineWidth',1, 'DisplayName','1:1 line');
    xlim([0 max_ax]); ylim([0 max_ax]);
    ax2 = gca;
    ax2.XTick = 0:2:ceil(max_ax);
    ax2.YTick = 0:2:ceil(max_ax);
    xlabel('IDW neighbourhood estimate  P_{IDW}  (mm)', 'FontWeight','normal', 'FontSize', label_fontsize);
    ylabel('Target station rain  P_{target}  (mm)', 'FontWeight','normal', 'FontSize', label_fontsize);
    set(gca, 'FontWeight','normal', 'FontSize', axis_fontsize);
    title(sprintf('Target vs neighbourhood IDW (ALL) - %s\nn = %d, flagged = %d', ...
        clabel, numel(idx_class), sum(is_flag)), 'FontWeight','normal', 'FontSize', title_fontsize);
    legend('Location','southeast', 'FontWeight','normal', 'FontSize', legend_fontsize);
    grid on; box on; axis square
    saveas(fig2, fullfile(plot_folder, sprintf('Scatter_%s_ALL.png', cname)));
    fig2.Visible = 'on';
    savefig(fig2, fullfile(plot_folder, sprintf('Scatter_%s_ALL.fig', cname)));
    close(fig2);

    % ================= اسکترهای جدا برای هر سطح کوانتایل (Q90/Q95/Q99) =================
    for qi = 1:numel(quantile_levels)
        p   = quantile_levels(qi);
        plabel = p*100;

        if one_sided
            idx_extreme = m_class > q_thr_up(qi);
        else
            idx_extreme = (m_class < q_thr_low(qi)) | (m_class > q_thr_up(qi));
        end

        if ~any(idx_extreme), continue; end

        xq = x_pidw(idx_extreme);
        yq = y_rain(idx_extreme);
        fq = is_flag(idx_extreme);

        fig3 = figure('Visible','off','Position',[100 100 800 750]);
        hold on
        scatter(xq, yq, 45, [0 0.35 0.75], 'filled', ...
            'MarkerFaceAlpha',0.75, 'DisplayName','Data');
        max_x = max(xq, [], 'omitnan') * 1.05;
        max_y = max(yq, [], 'omitnan') * 1.05;
        if isempty(max_x) || max_x <= 0 || isnan(max_x), max_x = 1; end
        if isempty(max_y) || max_y <= 0 || isnan(max_y), max_y = 1; end
        max_axq = max(max_x, max_y);
        plot([0 max_axq], [0 max_axq], 'k--', 'LineWidth',1, 'DisplayName','1:1 line');
        xlim([0 max_x]); ylim([0 max_y]);
        ax3 = gca;
        ax3.XTick = 0:2:ceil(max_x);
        ax3.YTick = 0:2:ceil(max_y);
        xlabel('IDW neighbourhood estimate  P_{IDW}  (mm)', 'FontWeight','normal', 'FontSize', label_fontsize);
        ylabel('Target station rain  P_{target}  (mm)', 'FontWeight','normal', 'FontSize', label_fontsize);
        set(gca, 'FontWeight','normal', 'FontSize', axis_fontsize);
        title(sprintf('Target vs neighbourhood IDW - %s - beyond Q%g\nn(extreme) = %d / %d total, flagged in this set = %d', ...
            clabel, plabel, sum(idx_extreme), numel(idx_class), sum(fq)), 'FontWeight','normal', 'FontSize', title_fontsize);
        legend('Location','southeast', 'FontWeight','normal', 'FontSize', legend_fontsize);
        grid on; box on
        saveas(fig3, fullfile(plot_folder, sprintf('Scatter_%s_Q%g.png', cname, plabel)));
        fig3.Visible = 'on';
        savefig(fig3, fullfile(plot_folder, sprintf('Scatter_%s_Q%g.fig', cname, plabel)));
        close(fig3);
    end

    % ================= NEW: اسکتر ترکیبی -- هر سه کوانتایل با رنگ‌های متفاوت در یک نمودار =================
    % به‌جای سه نمودار جدا (که بالاتر هم رسم شدند)، این یک نمودار اضافه
    % است که هر سه باند کوانتایل (Q90-Q95 / Q95-Q99 / بالاتر از Q99) را
    % هم‌زمان و با رنگ‌های متفاوت (مطابق همان رنگ‌های خطوط هیستوگرام/CDF
    % بالا) در یک محور واحد نشان می‌دهد.
    if one_sided
        band90 = m_class > q_thr_up(1) & m_class <= q_thr_up(2);   % بین Q90 و Q95
        band95 = m_class > q_thr_up(2) & m_class <= q_thr_up(3);   % بین Q95 و Q99
        band99 = m_class > q_thr_up(3);                             % بالاتر از Q99
    else
        band90 = (m_class < q_thr_low(1) & m_class >= q_thr_low(2)) | (m_class > q_thr_up(1) & m_class <= q_thr_up(2));
        band95 = (m_class < q_thr_low(2) & m_class >= q_thr_low(3)) | (m_class > q_thr_up(2) & m_class <= q_thr_up(3));
        band99 = (m_class < q_thr_low(3)) | (m_class > q_thr_up(3));
    end

    if any(band90) || any(band95) || any(band99)
        band_any = band90 | band95 | band99;

        fig4 = figure('Visible','off','Position',[100 100 800 750]);
        hold on
        scatter(x_pidw(band90), y_rain(band90), 40, quantile_colors(1,:), 'filled', ...
            'MarkerFaceAlpha',0.7, 'DisplayName','Q90-Q95');
        scatter(x_pidw(band95), y_rain(band95), 40, quantile_colors(2,:), 'filled', ...
            'MarkerFaceAlpha',0.7, 'DisplayName','Q95-Q99');
        scatter(x_pidw(band99), y_rain(band99), 55, quantile_colors(3,:), 'filled', ...
            'MarkerFaceAlpha',0.85, 'DisplayName','> Q99');

        max_x4 = max(x_pidw(band_any), [], 'omitnan') * 1.05;
        max_y4 = max(y_rain(band_any), [], 'omitnan') * 1.05;
        if isempty(max_x4) || max_x4 <= 0 || isnan(max_x4), max_x4 = 1; end
        if isempty(max_y4) || max_y4 <= 0 || isnan(max_y4), max_y4 = 1; end

        % ---- محدوده و گام تیک محور قائم -- فقط برای این نمودارهای جدید ----
        y_min4  = 0;
        yStep4  = 2;
        if cname == "Heavy"
            y_min4 = 10;     % برای بارش سنگین، محور قائم دقیقاً از 10 شروع شود
        elseif cname == "Moderate"
            y_min4 = 2.5;    % برای بارش متوسط، محور قائم دقیقاً از مرز پایین کلاس (2.5) شروع شود
            max_y4 = 10;     % سقف محور دقیقاً مرز بالای کلاس (بدون هیچ حاشیه‌ی اضافه)
        elseif cname == "Light"
            y_min4 = 0;      % مرز پایین کلاس Light
            max_y4 = 2.5;    % مرز بالای کلاس Light (بدون حاشیه‌ی اضافه)
            yStep4 = 0.5;    % تیک‌های نیم‌نیم: 0, 0.5, 1, 1.5, 2, 2.5
        end

        max_ax4 = max(max_x4, max_y4);   % بعد از اعمال override های بالا محاسبه می‌شود تا خط 1:1 هماهنگ بماند

        plot([0 max_ax4], [0 max_ax4], 'k--', 'LineWidth',1, 'DisplayName','1:1 line');
        xlim([0 max_x4]); ylim([y_min4 max_y4]);
        ax4 = gca;
        ax4.XTick = 0:2:ceil(max_x4);
        if cname == "Moderate"
            ax4.YTick = [2.5 4 5.5 7 8.5 10];   % تیک‌های ثابت و تمیز، بدون اعشار زشت
        else
            yticks_vec = y_min4:yStep4:max_y4;
            if isempty(yticks_vec) || abs(yticks_vec(end) - max_y4) > 1e-9
                yticks_vec = [yticks_vec, max_y4];
            end
            ax4.YTick = yticks_vec;
        end
        xlabel('IDW neighbourhood estimate  P_{IDW}  (mm)', 'FontWeight','normal', 'FontSize', label_fontsize);
        ylabel('Target station rain  P_{target}  (mm)', 'FontWeight','normal', 'FontSize', label_fontsize);
        set(gca, 'FontWeight','normal', 'FontSize', axis_fontsize);
        th4 = title(sprintf('Target vs neighbourhood IDW - %s - all quantile bands (Q90/Q95/Q99)\nn(Q90-Q95)=%d, n(Q95-Q99)=%d, n(>Q99)=%d', ...
            clabel, sum(band90), sum(band95), sum(band99)), 'FontWeight','normal', 'FontSize', title_fontsize);
        if ismember(cname, ["Heavy","Moderate","Light"])
            th4.Units = 'normalized';
            th4.Position(2) = th4.Position(2) + 0.015;   % کمی بالاتر بردن عنوان، دور از نقاط بالای نمودار
        end
        legend('Location','southeast', 'FontWeight','normal', 'FontSize', legend_fontsize);
        grid on; box on
        saveas(fig4, fullfile(plot_folder, sprintf('Scatter_%s_AllQuantiles.png', cname)));
        fig4.Visible = 'on';
        savefig(fig4, fullfile(plot_folder, sprintf('Scatter_%s_AllQuantiles.fig', cname)));
        close(fig4);
    end
end

fprintf('Diagnostic plots saved to: %s\n', plot_folder);

%% ---- 8. ساخت جدول‌های خروجی هر ایستگاه (ORIGINAL: Consorzio_hourly_Final_Table) ----
fprintf('Building final classified tables (ORIGINAL) with neighbour columns...\n');

Consorzio_hourly_Final_Table = table;
Consorzio_hourly_Final_Table.Station = Consorzio_hourly_Final(:,1);

Consorzio_hourly_Final_Table.Total          = cell(n_stations,1);
Consorzio_hourly_Final_Table.Artificial     = cell(n_stations,1);
Consorzio_hourly_Final_Table.FlagRepeat4h   = cell(n_stations,1);
Consorzio_hourly_Final_Table.Valid1         = cell(n_stations,1);
Consorzio_hourly_Final_Table.Flag200        = cell(n_stations,1);
Consorzio_hourly_Final_Table.Valid2         = cell(n_stations,1);
Consorzio_hourly_Final_Table.NotEvaluated   = cell(n_stations,1);
Consorzio_hourly_Final_Table.FlagNotAligned = cell(n_stations,1);
Consorzio_hourly_Final_Table.Valid3         = cell(n_stations,1);

%% ---- 8-NEW. ساخت جدول Q95: Consorzio_hourly_Final_Table_95 ----
fprintf('Building Consorzio_hourly_Final_Table_95 (Q95 rule applied to Light/Moderate/Heavy SARD)...\n');

Consorzio_hourly_Final_Table_95 = table;
Consorzio_hourly_Final_Table_95.Station = Consorzio_hourly_Final(:,1);

Consorzio_hourly_Final_Table_95.Total          = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.Artificial     = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.FlagRepeat4h   = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.Valid1         = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.Flag200        = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.Valid2         = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.NotEvaluated   = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.FlagNotAligned = cell(n_stations,1);
Consorzio_hourly_Final_Table_95.Valid3         = cell(n_stations,1);

for i = 1:n_stations
    d = Consorzio_hourly_Final{i,2};
    d.Artificial = artificial_matrix(:,i);
    d.NumValidNb = NumValidNb_mat(:,i);
    d.PIDW       = PIDW_mat(:,i);
    d.Class      = Class_mat(:,i);
    d.Metric     = Metric_Combined_mat(:,i);   % مقدار SARD (یا P_IDW برای NoRain) خود رکورد
    d.FlagAllNeighboursZero = FlagAllNeighboursZero_matrix(:,i);

    nb_col      = nb_col_all{i};
    nb_dist     = nb_dist_all{i};
    nb_xy_found = nb_xy_all{i};
    n_nb_f      = numel(nb_col);

    nb_names_found = strings(n_nb_f, 1);
    for k = 1:n_nb_f
        nb_names_found(k) = string(Consorzio_hourly_Final{nb_col(k), 1});
    end

    % ---------- نسخه‌ی ORIGINAL ----------
    d_total_o    = d;
    d_valid2_o   = d(Valid2_matrix(:,i), :);
    d_valid3_o   = d(Valid3_matrix(:,i), :);
    d_noteval_o  = d(NotEvaluated_matrix(:,i), :);
    d_notalign_o = d(FlagNotAligned_matrix(:,i), :);

    % ---------- نسخه‌ی Q95 ----------
    d_total_95    = d;
    d_valid2_95   = d(Valid2_matrix(:,i), :);          % Valid2 در هر دو نسخه یکسان است
    d_valid3_95   = d(Valid3_matrix_95(:,i), :);
    d_noteval_95  = d(NotEvaluated_matrix(:,i), :);    % NotEvaluated در هر دو نسخه یکسان است
    d_notalign_95 = d(FlagNotAligned_matrix_95(:,i), :);

    if n_nb_f > 0
        d_total_o    = addNeighbourColumns(d_total_o,    nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_valid2_o   = addNeighbourColumns(d_valid2_o,   nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_valid3_o   = addNeighbourColumns(d_valid3_o,   nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_noteval_o  = addNeighbourColumns(d_noteval_o,  nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_notalign_o = addNeighbourColumns(d_notalign_o, nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);

        d_total_95    = addNeighbourColumns(d_total_95,    nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_valid2_95   = addNeighbourColumns(d_valid2_95,   nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_valid3_95   = addNeighbourColumns(d_valid3_95,   nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_noteval_95  = addNeighbourColumns(d_noteval_95,  nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
        d_notalign_95 = addNeighbourColumns(d_notalign_95, nb_col, nb_dist, nb_names_found, nb_xy_found, rain_matrix, Valid2_matrix, time_centers, n_hours);
    end

    % ---- پر کردن جدول ORIGINAL ----
    Consorzio_hourly_Final_Table.Total{i}         = d_total_o;
    Consorzio_hourly_Final_Table.Artificial{i}    = d(artificial_matrix(:,i), :);
    Consorzio_hourly_Final_Table.FlagRepeat4h{i}  = d(FlagRepeat4h_matrix(:,i), :);
    Consorzio_hourly_Final_Table.Valid1{i}        = d(Valid1_matrix(:,i), :);
    Consorzio_hourly_Final_Table.Flag200{i}       = d(Flag200_matrix(:,i), :);
    Consorzio_hourly_Final_Table.Valid2{i}        = d_valid2_o;
    Consorzio_hourly_Final_Table.NotEvaluated{i}  = d_noteval_o;
    Consorzio_hourly_Final_Table.FlagNotAligned{i}= d_notalign_o;
    Consorzio_hourly_Final_Table.Valid3{i}        = d_valid3_o;

    % ---- پر کردن جدول Q95 ----
    Consorzio_hourly_Final_Table_95.Total{i}          = d_total_95;
    Consorzio_hourly_Final_Table_95.Artificial{i}     = d(artificial_matrix(:,i), :);
    Consorzio_hourly_Final_Table_95.FlagRepeat4h{i}   = d(FlagRepeat4h_matrix(:,i), :);
    Consorzio_hourly_Final_Table_95.Valid1{i}         = d(Valid1_matrix(:,i), :);
    Consorzio_hourly_Final_Table_95.Flag200{i}        = d(Flag200_matrix(:,i), :);
    Consorzio_hourly_Final_Table_95.Valid2{i}         = d_valid2_95;
    Consorzio_hourly_Final_Table_95.NotEvaluated{i}   = d_noteval_95;
    Consorzio_hourly_Final_Table_95.FlagNotAligned{i} = d_notalign_95;
    Consorzio_hourly_Final_Table_95.Valid3{i}         = d_valid3_95;
end

%% ---- 8b. Cumulative annual rainfall per station: Total vs Valid3 (هر دو نسخه) ----
fprintf('Computing cumulative annual rainfall (Total vs Valid3), both versions...\n');

Consorzio_hourly_Final_Table.CumulativeTotal_mm  = zeros(n_stations,1);
Consorzio_hourly_Final_Table.CumulativeValid3_mm = zeros(n_stations,1);
Consorzio_hourly_Final_Table.RemovedByQC_mm      = zeros(n_stations,1);
Consorzio_hourly_Final_Table.RemovedByQC_pct     = zeros(n_stations,1);

Consorzio_hourly_Final_Table_95.CumulativeTotal_mm  = zeros(n_stations,1);
Consorzio_hourly_Final_Table_95.CumulativeValid3_mm = zeros(n_stations,1);
Consorzio_hourly_Final_Table_95.RemovedByQC_mm      = zeros(n_stations,1);
Consorzio_hourly_Final_Table_95.RemovedByQC_pct     = zeros(n_stations,1);

for i = 1:n_stations
    cum_total     = sum(Consorzio_hourly_Final_Table.Total{i}.rain,  'omitnan');
    cum_valid3_o  = sum(Consorzio_hourly_Final_Table.Valid3{i}.rain, 'omitnan');
    cum_valid3_95 = sum(Consorzio_hourly_Final_Table_95.Valid3{i}.rain, 'omitnan');

    Consorzio_hourly_Final_Table.CumulativeTotal_mm(i)  = cum_total;
    Consorzio_hourly_Final_Table.CumulativeValid3_mm(i) = cum_valid3_o;
    Consorzio_hourly_Final_Table.RemovedByQC_mm(i)      = cum_total - cum_valid3_o;

    Consorzio_hourly_Final_Table_95.CumulativeTotal_mm(i)  = cum_total;
    Consorzio_hourly_Final_Table_95.CumulativeValid3_mm(i) = cum_valid3_95;
    Consorzio_hourly_Final_Table_95.RemovedByQC_mm(i)      = cum_total - cum_valid3_95;

    if cum_total > 0
        Consorzio_hourly_Final_Table.RemovedByQC_pct(i)    = 100 * (cum_total - cum_valid3_o)  / cum_total;
        Consorzio_hourly_Final_Table_95.RemovedByQC_pct(i) = 100 * (cum_total - cum_valid3_95) / cum_total;
    else
        Consorzio_hourly_Final_Table.RemovedByQC_pct(i)    = 0;
        Consorzio_hourly_Final_Table_95.RemovedByQC_pct(i) = 0;
    end
end

%% ---- 9. گزارش خلاصه (ORIGINAL) ----
fprintf('\n========== Quality Control Summary (ORIGINAL, L/M/H not SARD-flagged) ==========\n');
fprintf('%-30s %6s %6s %8s %7s %7s %7s %11s %11s %7s\n', ...
    'Station','Total','Artif','Rep4h','V1','>200','V2','NotEval','NotAlign','V3');
fprintf('%s\n', repmat('-',1,105));
for i = 1:n_stations
    sn = char(Consorzio_hourly_Final_Table.Station{i});
    if numel(sn)>29, sn=sn(1:29); end
    fprintf('%-30s %6d %6d %8d %7d %7d %7d %11d %11d %7d\n', sn, ...
        height(Consorzio_hourly_Final_Table.Total{i}), ...
        height(Consorzio_hourly_Final_Table.Artificial{i}), ...
        height(Consorzio_hourly_Final_Table.FlagRepeat4h{i}), ...
        height(Consorzio_hourly_Final_Table.Valid1{i}), ...
        height(Consorzio_hourly_Final_Table.Flag200{i}), ...
        height(Consorzio_hourly_Final_Table.Valid2{i}), ...
        height(Consorzio_hourly_Final_Table.NotEvaluated{i}), ...
        height(Consorzio_hourly_Final_Table.FlagNotAligned{i}), ...
        height(Consorzio_hourly_Final_Table.Valid3{i}));
end

fprintf('\n========== Quality Control Summary (Q95, L/M/H flagged if SARD >= class Q95) ==========\n');
fprintf('%-30s %6s %6s %8s %7s %7s %7s %11s %11s %7s\n', ...
    'Station','Total','Artif','Rep4h','V1','>200','V2','NotEval','NotAlign','V3');
fprintf('%s\n', repmat('-',1,105));
for i = 1:n_stations
    sn = char(Consorzio_hourly_Final_Table_95.Station{i});
    if numel(sn)>29, sn=sn(1:29); end
    fprintf('%-30s %6d %6d %8d %7d %7d %7d %11d %11d %7d\n', sn, ...
        height(Consorzio_hourly_Final_Table_95.Total{i}), ...
        height(Consorzio_hourly_Final_Table_95.Artificial{i}), ...
        height(Consorzio_hourly_Final_Table_95.FlagRepeat4h{i}), ...
        height(Consorzio_hourly_Final_Table_95.Valid1{i}), ...
        height(Consorzio_hourly_Final_Table_95.Flag200{i}), ...
        height(Consorzio_hourly_Final_Table_95.Valid2{i}), ...
        height(Consorzio_hourly_Final_Table_95.NotEvaluated{i}), ...
        height(Consorzio_hourly_Final_Table_95.FlagNotAligned{i}), ...
        height(Consorzio_hourly_Final_Table_95.Valid3{i}));
end

fprintf('\n========== Cumulative Annual Rainfall (ORIGINAL: Total vs Valid3) ==========\n');
fprintf('%-30s %14s %14s %12s %9s\n', 'Station','Total (mm)','Valid3 (mm)','Removed(mm)','Removed(%)');
fprintf('%s\n', repmat('-',1,85));
for i = 1:n_stations
    sn = char(Consorzio_hourly_Final_Table.Station{i});
    if numel(sn)>29, sn=sn(1:29); end
    fprintf('%-30s %14.2f %14.2f %12.2f %8.2f%%\n', sn, ...
        Consorzio_hourly_Final_Table.CumulativeTotal_mm(i), ...
        Consorzio_hourly_Final_Table.CumulativeValid3_mm(i), ...
        Consorzio_hourly_Final_Table.RemovedByQC_mm(i), ...
        Consorzio_hourly_Final_Table.RemovedByQC_pct(i));
end

fprintf('\n========== Cumulative Annual Rainfall (Q95: Total vs Valid3) ==========\n');
fprintf('%-30s %14s %14s %12s %9s\n', 'Station','Total (mm)','Valid3 (mm)','Removed(mm)','Removed(%)');
fprintf('%s\n', repmat('-',1,85));
for i = 1:n_stations
    sn = char(Consorzio_hourly_Final_Table_95.Station{i});
    if numel(sn)>29, sn=sn(1:29); end
    fprintf('%-30s %14.2f %14.2f %12.2f %8.2f%%\n', sn, ...
        Consorzio_hourly_Final_Table_95.CumulativeTotal_mm(i), ...
        Consorzio_hourly_Final_Table_95.CumulativeValid3_mm(i), ...
        Consorzio_hourly_Final_Table_95.RemovedByQC_mm(i), ...
        Consorzio_hourly_Final_Table_95.RemovedByQC_pct(i));
end

%% ---- 10. ذخیره ----
output_folder = 'Output_Consorzio_hourly_v3';
if ~exist(output_folder,'dir'), mkdir(output_folder); end

thresholds_used = struct( ...
    'Hourly_threshold',          Hourly_threshold, ...
    'Repeat_run_length',         Repeat_run_length, ...
    'Min_valid_neighbours',      Min_valid_neighbours, ...
    'Light_upper_mm',            Light_upper_mm, ...
    'Moderate_upper_mm',         Moderate_upper_mm, ...
    'NoRain_fixed_threshold_mm', NoRain_fixed_threshold_mm, ...
    'AllNeighboursZero_target_threshold_mm', AllNeighboursZero_target_threshold_mm, ...
    'thr_norain_upper',          thr_norain_upper);

save(fullfile(output_folder,'Consorzio_hourly_Final_Table_v2.mat'), ...
     'Consorzio_hourly_Final_Table', 'thresholds_used');

thresholds_used_95 = struct( ...
    'Hourly_threshold',          Hourly_threshold, ...
    'Repeat_run_length',         Repeat_run_length, ...
    'Min_valid_neighbours',      Min_valid_neighbours, ...
    'Light_upper_mm',            Light_upper_mm, ...
    'Moderate_upper_mm',         Moderate_upper_mm, ...
    'NoRain_fixed_threshold_mm', NoRain_fixed_threshold_mm, ...
    'AllNeighboursZero_target_threshold_mm', AllNeighboursZero_target_threshold_mm, ...
    'thr_norain_upper',          thr_norain_upper, ...
    'Q95_level',                 Q95_level, ...
    'thr_light_q95',             thr_light_q95, ...
    'thr_moderate_q95',          thr_moderate_q95, ...
    'thr_heavy_q95',             thr_heavy_q95);

save(fullfile(output_folder,'Consorzio_hourly_Final_Table_95.mat'), ...
     'Consorzio_hourly_Final_Table_95', 'thresholds_used_95');

fprintf('\nSaved: Consorzio_hourly_Final_Table_v2.mat  (ORIGINAL, L/M/H not SARD-flagged)\n');
fprintf('Saved: Consorzio_hourly_Final_Table_95.mat  (Q95 rule applied to L/M/H SARD)\n');

output_stations_folder = fullfile(output_folder, 'PerStation');
if ~exist(output_stations_folder,'dir'), mkdir(output_stations_folder); end

output_stations_folder_95 = fullfile(output_folder, 'PerStation_Q95');
if ~exist(output_stations_folder_95,'dir'), mkdir(output_stations_folder_95); end

for i = 1:n_stations
    raw_name  = char(Consorzio_hourly_Final_Table.Station{i});
    safe_name = regexprep(raw_name, '[^\w]', '_');
    safe_name = regexprep(safe_name, '_+', '_');
    safe_name = strtrim(safe_name);

    % ---- ORIGINAL ----
    Station_Total                = Consorzio_hourly_Final_Table.Total{i};
    Station_Artificial            = Consorzio_hourly_Final_Table.Artificial{i};
    Station_FlagRepeat4h          = Consorzio_hourly_Final_Table.FlagRepeat4h{i};
    Station_Valid1                = Consorzio_hourly_Final_Table.Valid1{i};
    Station_Flag200                = Consorzio_hourly_Final_Table.Flag200{i};
    Station_Valid2                  = Consorzio_hourly_Final_Table.Valid2{i};
    Station_NotEvaluated            = Consorzio_hourly_Final_Table.NotEvaluated{i};
    Station_FlagNotAligned          = Consorzio_hourly_Final_Table.FlagNotAligned{i};
    Station_Valid3                  = Consorzio_hourly_Final_Table.Valid3{i};
    Station_CumulativeTotal_mm      = Consorzio_hourly_Final_Table.CumulativeTotal_mm(i);
    Station_CumulativeValid3_mm     = Consorzio_hourly_Final_Table.CumulativeValid3_mm(i);
    Station_RemovedByQC_mm          = Consorzio_hourly_Final_Table.RemovedByQC_mm(i);
    Station_RemovedByQC_pct         = Consorzio_hourly_Final_Table.RemovedByQC_pct(i);

    save_path = fullfile(output_stations_folder, [safe_name '.mat']);
    save(save_path, ...
        'Station_Total','Station_Artificial','Station_FlagRepeat4h', ...
        'Station_Valid1','Station_Flag200','Station_Valid2', ...
        'Station_NotEvaluated','Station_FlagNotAligned','Station_Valid3', ...
        'Station_CumulativeTotal_mm','Station_CumulativeValid3_mm', ...
        'Station_RemovedByQC_mm','Station_RemovedByQC_pct');

    % ---- Q95 ----
    Station_Total_95                = Consorzio_hourly_Final_Table_95.Total{i};
    Station_Artificial_95           = Consorzio_hourly_Final_Table_95.Artificial{i};
    Station_FlagRepeat4h_95         = Consorzio_hourly_Final_Table_95.FlagRepeat4h{i};
    Station_Valid1_95               = Consorzio_hourly_Final_Table_95.Valid1{i};
    Station_Flag200_95              = Consorzio_hourly_Final_Table_95.Flag200{i};
    Station_Valid2_95               = Consorzio_hourly_Final_Table_95.Valid2{i};
    Station_NotEvaluated_95         = Consorzio_hourly_Final_Table_95.NotEvaluated{i};
    Station_FlagNotAligned_95       = Consorzio_hourly_Final_Table_95.FlagNotAligned{i};
    Station_Valid3_95               = Consorzio_hourly_Final_Table_95.Valid3{i};
    Station_CumulativeTotal_mm_95   = Consorzio_hourly_Final_Table_95.CumulativeTotal_mm(i);
    Station_CumulativeValid3_mm_95  = Consorzio_hourly_Final_Table_95.CumulativeValid3_mm(i);
    Station_RemovedByQC_mm_95       = Consorzio_hourly_Final_Table_95.RemovedByQC_mm(i);
    Station_RemovedByQC_pct_95      = Consorzio_hourly_Final_Table_95.RemovedByQC_pct(i);

    save_path_95 = fullfile(output_stations_folder_95, [safe_name '.mat']);
    save(save_path_95, ...
        'Station_Total_95','Station_Artificial_95','Station_FlagRepeat4h_95', ...
        'Station_Valid1_95','Station_Flag200_95','Station_Valid2_95', ...
        'Station_NotEvaluated_95','Station_FlagNotAligned_95','Station_Valid3_95', ...
        'Station_CumulativeTotal_mm_95','Station_CumulativeValid3_mm_95', ...
        'Station_RemovedByQC_mm_95','Station_RemovedByQC_pct_95');

    fprintf('Saved station %d/%d: %s (ORIGINAL + Q95)\n', i, n_stations, safe_name);
end

fprintf('\nAll files saved to: %s\n', output_folder);
fprintf('Diagnostic plots saved to: %s\n', plot_folder);
toc

%% ============================================================
%  LOCAL FUNCTIONS
%% ============================================================

%% ============================================================
%  Local functions (moved to the end of the script so it also runs
%  on MATLAB releases older than R2024a)
%% ============================================================
function T_out = addNeighbourColumns(T_sub, nb_col, nb_dist, ...
                                      nb_names_found, nb_xy, ...
                                      rain_matrix, valid2_matrix, ...
                                      time_centers, n_hours)
    n_nb  = numel(nb_col);
    n_row = height(T_sub);
    [~, time_idx] = ismember(T_sub.time, time_centers);
    T_out = T_sub;
    for k = 1:n_nb
        col_idx        = nb_col(k);
        nb_rain_sel    = rain_matrix(time_idx, col_idx);
        nb_valid2_sel  = valid2_matrix(time_idx, col_idx);
        prefix = sprintf('Nb%d_', k);
        T_out.(sprintf('%sName',    prefix)) = repmat(string(nb_names_found(k)), n_row, 1);
        T_out.(sprintf('%sDistance',prefix)) = repmat(nb_dist(k), n_row, 1);
        T_out.(sprintf('%sX',       prefix)) = repmat(nb_xy(k,1), n_row, 1);
        T_out.(sprintf('%sY',       prefix)) = repmat(nb_xy(k,2), n_row, 1);
        T_out.(sprintf('%sRain',    prefix)) = nb_rain_sel;
        T_out.(sprintf('%sValid2',  prefix)) = nb_valid2_sel;
    end
end
