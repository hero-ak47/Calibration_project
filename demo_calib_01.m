%% ==============================================================
%  DEMO: Thuat toan Calibration pha cho tuyen phat Massive MIMO 5G
%  Tham khao: 
%   - Tran et al., "Real-time Calibration for Digital Beamforming in 
%     5G Systems With Experiments on Testbed", IEEE ATC 2022 
%     (DOI: 10.1109/ATC55345.2022.9943013)
%   - Bang doc quyen sang che VN so 1-0052769: "Phuong phap dieu chinh 
%     pha tin hieu cho tuyen phat ung dung cho he thong thu phat song 
%     vo tuyen 5G"
%
%  MUC TIEU DEMO
%  --------------
%  Ca bai bao va bang sang che deu mo ta thuat toan qua 3 buoc:
%    B1: Phat - thu tin hieu tham chieu (pilot ZC) qua coupler/duong 
%        phan hoi noi bo.
%    B2: Uoc luong dap ung kenh H tren tung sub-carrier, tu do suy ra 
%        do lech pha gay boi phan cung.
%    B3: Tinh he so bu pha va ap dung cho du lieu phat.
%
%  Trong thuc te, do lech pha do phan cung gay ra gom 2 thanh phan:
%    (a) GROUP DELAY  - dich pha TUYEN TINH theo sub-carrier k (do do 
%        dai duong truyen/cap khac nhau giua cac chain)
%    (b) FINE PHASE   - dich pha HANG SO, khong phu thuoc k (do lech 
%        pha LO/mixer/PA)
%  Neu chi bu MOT thanh phan (vi du chi do pha tai 1 sub-carrier trung 
%  tam roi ap dung cho ca bang thong), phan GROUP DELAY se khong duoc 
%  bu -> loi con lai tang dan khi ra xa sub-carrier do -> beam squint.
%  Day chinh la buoc "an" ma nguoi dung nghi la tac gia co lam nhung 
%  khong noi ro trong phan mo ta rut gon.
%
%  Demo nay:
%   1. Mo phong N_ant chuoi RF voi loi phan cung "that" (tau_n, phi0_n)
%   2. Uoc luong 2 buoc: Group delay (tho) truoc, Fine phase (tinh) sau
%      - dung dung cong thuc (2)-(5) trong bai bao IEEE ATC 2022
%   3. So sanh 3 kich ban bu pha: Khong calib / Chi Fine-phase / Day du
%   4. Ve sai pha con lai theo sub-carrier va bup song beamforming
% ================================================================

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

fprintf('--- GROUND TRUTH (loi phan cung that, chi de doi chieu) ---\n');
for a = 1:N_ant
    fprintf('  Chain %d: tau = %6.2f mau | phi0 = %7.2f do\n', ...
        a, true_tau_n(a), true_phi0_n(a)*180/pi);
end

%% 4. BUOC 1 (mo phong): PHAT-THU TIN HIEU THAM CHIEU QUA COUPLER, L LAN LAP
k_idx = (0:N_sub-1);              % chi so sub-carrier (0-based, giong cong thuc paper)
Ybar  = zeros(N_ant, N_sub);      % Y_n^k trung binh sau L lan (cong thuc (3))

for a = 1:N_ant
    acc = zeros(1, N_sub);
    H_true = exp(1j*true_phi0_n(a)) .* exp(1j*2*pi*k_idx*true_tau_n(a)/M);
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

for a = 1:N_ant
    h = Hhat(a,:);
    % --- (a) UOC LUONG GROUP DELAY (THO) - dung cong thuc (5) trong paper ---
    prod_term  = h(1:end-delta_k) .* conj(h(1+delta_k:end));
    phase_diff = angle(mean(prod_term));
    tau_hat(a) = -M/(2*pi*delta_k) * phase_diff;

    % --- (b) BU GROUP DELAY, UOC LUONG PHA DU (FINE) ---
    h_comp    = h .* exp(-1j*2*pi*k_idx*tau_hat(a)/M);
    phi0_hat(a) = angle(mean(h_comp));
end

fprintf('\n--- KET QUA UOC LUONG (so sanh voi ground truth) ---\n');
for a = 1:N_ant
    fprintf('  Chain %d: tau_hat = %6.2f mau (loi %+.3f) | phi0_hat = %7.2f do (loi %+.2f do)\n', ...
        a, tau_hat(a), tau_hat(a)-true_tau_n(a), ...
        phi0_hat(a)*180/pi, (phi0_hat(a)-true_phi0_n(a))*180/pi);
end

%% 6. BUOC 3: TINH 3 BO TRONG SO CALIBRATION DE SO SANH
% (1) Khong calib
W_none = ones(N_ant, N_sub);

% (2) CHI bu Fine-phase tai 1 sub-carrier trung tam, BO QUA group delay
%     -> mo phong dung "buoc bi an" ma nguoi dung nghi tac gia co lam
k0          = round(N_sub/2);
phi_center  = angle(Hhat(:, k0));              % N_ant x 1
W_fineOnly  = repmat(exp(-1j*phi_center), 1, N_sub);

% (3) Bu DAY DU: group delay + fine phase (dung phuong phap paper/patent)
W_full = exp(-1j*(2*pi*k_idx.*tau_hat.'/M + phi0_hat.'));

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
        H_true_k = exp(1j*true_phi0_n.') .* exp(1j*2*pi*(k-1)*true_tau_n.'/M);
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

fprintf('\nHoan tat demo. Xem 2 figure: (1) sai pha con lai theo sub-carrier, ');
fprintf('(2) bup song phat tai 3 vi tri trong bang thong.\n');
