%% Calibration bu pha va bien do
% PHASE ( lech thoi gian tau)
% - Gom 2 thanh phan Group delay va Fine delay
% --Group delay: + Duoc tinh thong qua trung binh do lech pha cua cac dap ung kenh cua tung subcarrier
%                + Voi cung 1 do lech thoi gian tau, do lech pha tren tung song mang con ( Yk=Xk⋅e^(−j2πfk​τ)) tuyen tinh theo k. Nghia la cac song mang con
%                  co tan so khac nhau nen do lech pha se ti le voi k (tan so cang cao thi lech pha cang lon)
%--Fine delay : + Day la do lech pha sau khi da bu tho (van con phan group delay chua bu chinh xac)
%               + Duoc tinh thong qua goc pha trung binh cua tin hieu phat sau khi da bu group delay. Cac tin hieu sau khi bu group delay da nam gan voi tin hieu ly
%                 tuong, viec tinh trung binh cong se lam do lech bu tru nhau, cang sat toi tin hieu ly tuong hon.

% AMPLITUDE : Duoc tinh thong qua trung binh do lech bien do cua cac dap ung kenh cua tung subcarrier
clear; clc; close all;
rng(7);   % seed co dinh de ket qua lap lai duoc

%% 1. THAM SO HE THONG (rut gon tu Bang 2 trong sang che & Section III paper)
N_ant    = 8;          % so luong kenh anten phat (RF chain), giong RRU 8T8R
N_sub    = 256;        % so sub-carrier dung cho calib (rut gon tu 3276 de demo nhanh)
M        = 256;        % kich thuoc IFFT/FFT
L        = 6;          % so lan lap phat-thu de trung binh giam nhieu (giong Bang 2)
SNR_dB   = 25;         % SNR tin hieu tham chieu tren duong phan hoi coupler
Ncell_ID = 5;          % ID de sinh chuoi ZC (demo, giong tham so trong sang che)

d_over_lambda    = 0.5;    % khoang cach anten = lambda/2
theta_target_deg = -30;    % goc bup song mong muon (giong Hinh 12 trong sang che)

%% 2. SINH CHUOI PILOT THAM CHIEU THEO ZADOFF-CHU (Buoc 1)
u = 25;                                  % root index cua chuoi ZC
n = (0:N_sub-1);
zc_base = exp(-1j*pi*u*n.*(n+1)/N_sub);  % chuoi ZC co ban (CAZAC)
X = exp(1j*Ncell_ID*n) .* zc_base;       % X_n^k: tin hieu tham chieu mien tan so
X = X ./ abs(X);                         % chuan hoa bien do = 1

%% 3. MO HINH LOI PHAN CUNG "THAT" CUA TUNG CHUOI RF (ground truth, an voi
%     thuat toan uoc luong - chi dung de kiem chung ket qua)
true_tau_n  = [0, 1.7, -2.3, 3.1, 0.6, -1.2, 2.8, -0.4];        % mau (group delay)
true_phi0_n = [0, 35, -60, 120, -150, 80, -25, 170] * pi/180;    % rad (fine phase)
% LECH BIEN DO: do dung sai gain PA, suy hao cap/connector khac nhau giua
% cac chain - la mot hang so (khong phu thuoc k) rieng cho tung anten,
% tuong tu ban chat voi fine phase nhung tac dong len BIEN DO thay vi PHA
true_A_n    = [1.0, 0.82, 1.18, 0.90, 1.10, 0.78, 1.22, 0.95];  % (khong don vi, ty le so voi chain 1)

fprintf('--- GROUND TRUTH (loi phan cung that, chi de doi chieu) ---\n');
for a = 1:N_ant
    fprintf('  Chain %d: tau = %6.2f mau | phi0 = %7.2f do | A = %5.2f\n', ...
        a, true_tau_n(a), true_phi0_n(a)*180/pi, true_A_n(a));
end

%% 4. BUOC 1 (mo phong): PHAT-THU TIN HIEU THAM CHIEU QUA COUPLER, L LAN LAP
k_idx = (0:N_sub-1);              % chi so sub-carrier (0-based, giong cong thuc paper)
Ybar  = zeros(N_ant, N_sub);      % Y_n^k trung binh sau L lan (cong thuc (3))

for a = 1:N_ant
    acc = zeros(1, N_sub);
    % Dap ung kenh THAT gio gom ca 3 thanh phan: bien do A_n, fine phase
    % phi0_n (hang so), va group delay tau_n (tuyen tinh theo k)
    H_true = true_A_n(a) * exp(1j*true_phi0_n(a)) .* exp(1j*2*pi*k_idx*true_tau_n(a)/M);
    for it = 1:L
        noise = (randn(1,N_sub) + 1j*randn(1,N_sub))/sqrt(2) * 10^(-SNR_dB/20);
        Y = X .* H_true + noise;
        acc = acc + Y;
    end
    Ybar(a,:) = acc / L;
end

%% 5. BUOC 2: UOC LUONG H, TACH GROUP DELAY (THO) & FINE PHASE (TINH)
Hhat = Ybar ./ X;      % H_n^k = Ybar / X   (cong thuc (2))

delta_k  = 8;                      % khoang cach sub-carrier de uoc luong group delay
tau_hat  = zeros(1, N_ant);
phi0_hat = zeros(1, N_ant);
A_hat    = zeros(1, N_ant);        % <-- MOI: uoc luong bien do tung chain

for a = 1:N_ant
    h = Hhat(a,:);
    % --- (a) UOC LUONG GROUP DELAY (THO) - dung cong thuc (5) trong paper ---
    prod_term  = h(1:end-delta_k) .* conj(h(1+delta_k:end));
    phase_diff = angle(mean(prod_term));
    tau_hat(a) = -M/(2*pi*delta_k) * phase_diff;

    % --- (b) BU GROUP DELAY, UOC LUONG PHA DU (FINE) VA BIEN DO ---
    h_comp    = h .* exp(-1j*2*pi*k_idx*tau_hat(a)/M);
    % Sau khi bu group delay, h_comp[k] ~ A_n * e^{j*phi0_n} (hang so
    % voi moi k, chi con lech boi nhieu) -> trung binh phasor (coherent
    % averaging) cho ca goc (angle) LAN do lon (abs) cung mot cong thuc:
    phi0_hat(a) = angle(mean(h_comp));
    A_hat(a)    = abs(mean(h_comp));      % <-- MOI: uoc luong bien do
end

fprintf('\n--- KET QUA UOC LUONG (so sanh voi ground truth) ---\n');
for a = 1:N_ant
    fprintf('  Chain %d: tau_hat=%6.2f mau (loi %+.3f) | phi0_hat=%7.2f do (loi %+.2f do) | A_hat=%5.3f (loi %+.2f%%)\n', ...
        a, tau_hat(a), tau_hat(a)-true_tau_n(a), ...
        phi0_hat(a)*180/pi, (phi0_hat(a)-true_phi0_n(a))*180/pi, ...
        A_hat(a), 100*(A_hat(a)-true_A_n(a))/true_A_n(a));
end

%% 6. BUOC 3: TINH 3 BO TRONG SO CALIBRATION DE SO SANH
% (1) Khong calib
W_none = ones(N_ant, N_sub);

% (2) CHI bu Fine-phase tai 1 sub-carrier trung tam, BO QUA group delay
%     -> mo phong dung "buoc bi an" ma nguoi dung nghi tac gia co lam
k0          = round(N_sub/2);
phi_center  = angle(Hhat(:, k0));              % N_ant x 1
W_fineOnly  = repmat(exp(-1j*phi_center), 1, N_sub);

% (3a) Bu DAY DU nhung CHI PHA: group delay + fine phase (giong paper/patent
%      goc, CHUA bu bien do) - dung de so sanh rieng anh huong cua bien do
W_full_phaseOnly = exp(-1j*(2*pi*k_idx.*tau_hat.'/M + phi0_hat.'));

% (3b) Bu DAY DU CA PHA VA BIEN DO (nang cap so voi paper/patent goc):
%      nhan them he so 1/A_hat de "keo" moi chain ve cung mot bien do
%      chuan hoa - tan dung luon |Hhat| von da co san, khong can do them
W_full = (1 ./ A_hat.') .* W_full_phaseOnly;

%% 7. SAI PHA CON LAI THEO SUB-CARRIER SAU KHI AP DUNG TUNG W
resid_none     = angle(Hhat .* W_none)     * 180/pi;
resid_fineOnly = angle(Hhat .* W_fineOnly) * 180/pi;
resid_full     = angle(Hhat .* W_full)     * 180/pi;

figure('Name','Sai pha con lai theo sub-carrier');
subplot(3,1,1); plot(k_idx, resid_none.');     ylim([-200 200]); grid on;
title('Khong calib'); ylabel('Pha (do)');
subplot(3,1,2); plot(k_idx, resid_fineOnly.'); ylim([-200 200]); grid on;
title('Chi bu Fine-Phase tai 1 sub-carrier (BO QUA Group Delay)'); ylabel('Pha (do)');
subplot(3,1,3); plot(k_idx, resid_full.');     ylim([-200 200]); grid on;
title('Bu day du: Group Delay + Fine Phase'); ylabel('Pha (do)'); xlabel('Chi so sub-carrier k');

rms_none     = sqrt(mean(resid_none(:).^2));
rms_fineOnly = sqrt(mean(resid_fineOnly(:).^2));
rms_full     = sqrt(mean(resid_full(:).^2));
fprintf('\n--- RMS SAI PHA TREN TOAN BANG THONG ---\n');
fprintf('  Khong calib          : %.2f do\n', rms_none);
fprintf('  Chi Fine-phase       : %.2f do\n', rms_fineOnly);
fprintf('  Group delay + Fine   : %.2f do\n', rms_full);

%% 8. MO PHONG BUP SONG (ARRAY FACTOR) TAI 3 VI TRI SUB-CARRIER KHAC NHAU
theta_target  = theta_target_deg*pi/180;
n_ant_idx     = (0:N_ant-1).';
a_theta       = exp(-1j*2*pi*d_over_lambda*n_ant_idx*sin(theta_target));  % steering vector

theta_scan_deg = -90:0.2:90;
theta_scan     = theta_scan_deg*pi/180;

k_eval_list = [1, round(N_sub/2), N_sub];   % dau bang - giua bang - cuoi bang
titles = {'Dau bang thong (k=1)','Giua bang thong','Cuoi bang thong (k=N)'};

figure('Name','Bup song phat theo goc quet - so sanh 3 kich ban calib');
for ik = 1:numel(k_eval_list)
    k = k_eval_list(ik);
    subplot(1,3,ik); hold on;
    for caseIdx = 1:3
        switch caseIdx
            case 1, Wk = W_none(:,k);     lbl = 'Khong calib';
            case 2, Wk = W_fineOnly(:,k); lbl = 'Chi Fine-phase';
            case 3, Wk = W_full(:,k);     lbl = 'Group delay + Fine';
        end
        % Dap ung phan cung THAT tai sub-carrier k (dung de kiem chung)
        % - gio da gom ca true_A_n (bien do that) de bup song phan anh
        %   dung ca 3 loai loi phan cung: delay + fine phase + bien do
        H_true_k = true_A_n.' .* exp(1j*true_phi0_n.') .* exp(1j*2*pi*(k-1)*true_tau_n.'/M);
        g = H_true_k .* Wk .* a_theta;         % dap ung tong hop moi anten
        AFpat = zeros(size(theta_scan));
        for it = 1:numel(theta_scan)
            sv = exp(-1j*2*pi*d_over_lambda*n_ant_idx*sin(theta_scan(it)));
            AFpat(it) = abs(sum(g .* conj(sv)))^2;
        end
        AFpat_dB = 10*log10(AFpat/max(AFpat) + eps);
        plot(theta_scan_deg, AFpat_dB, 'DisplayName', lbl, 'LineWidth', 1.3);
    end
    ylim([-40 0]); grid on; xlabel('Goc (do)'); ylabel('Cong suat chuan hoa (dB)');
    title(titles{ik}); legend('Location','south');
end
main_title = sprintf('Bup song huong %.0f do - so sanh cac phuong phap calib', theta_target_deg);
if exist('sgtitle', 'builtin') || exist('sgtitle', 'file')
    sgtitle(main_title);   % MATLAB R2018b+
else
    annotation('textbox', [0 0.94 1 0.05], 'String', main_title, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

%% 9. RIENG BIET: ANH HUONG CUA VIEC BU/KHONG BU BIEN DO LEN NULL & SIDELOBE
% O day CO DINH pha da bu day du (group delay + fine phase, dung nhu nhau
% o ca 2 truong hop) - CHI thay doi co bu bien do hay khong - de tach
% rieng dung anh huong cua bien do, khong lan voi anh huong cua pha.
k_center = round(N_sub/2);
H_true_center = true_A_n.' .* exp(1j*true_phi0_n.') .* exp(1j*2*pi*(k_center-1)*true_tau_n.'/M);

figure('Name','Anh huong cua bu bien do len sidelobe/null (da bu pha day du)');
hold on;
g_noAmp   = H_true_center .* W_full_phaseOnly(:,k_center) .* a_theta;   % pha dung, bien do CHUA bu
g_withAmp = H_true_center .* W_full(:,k_center)           .* a_theta;   % pha + bien do deu da bu
labels9 = {'Da bu pha, CHUA bu bien do','Da bu ca pha VA bien do'};
gs = {g_noAmp, g_withAmp};
for c = 1:2
    AFpat = zeros(size(theta_scan));
    for it = 1:numel(theta_scan)
        sv = exp(-1j*2*pi*d_over_lambda*n_ant_idx*sin(theta_scan(it)));
        AFpat(it) = abs(sum(gs{c} .* conj(sv)))^2;
    end
    plot(theta_scan_deg, 10*log10(AFpat/max(AFpat)+eps), 'LineWidth', 1.3, 'DisplayName', labels9{c});
end
ylim([-40 0]); grid on; legend('Location','south');
xlabel('Goc (do)'); ylabel('Cong suat chuan hoa (dB)');
title('Bu bien do khong doi HUONG bup song, chi lam null sau hon / sidelobe deu hon');

fprintf('\n--- SO SANH RIENG ANH HUONG CUA BU BIEN DO (pha da bu day du ca 2 truong hop) ---\n');
fprintf('  (Xem figure 3: huong dinh bup song giu nguyen -%d do o ca 2, chi khac null/sidelobe)\n', abs(theta_target_deg));

fprintf('\nHoan tat demo. Xem 3 figure: (1) sai pha con lai theo sub-carrier, ');
fprintf('(2) bup song phat tai 3 vi tri trong bang thong (da gom ca loi bien do), ');
fprintf('(3) anh huong rieng cua viec bu bien do len null/sidelobe.\n');
