%% Healthy vs. elastic-aneurysm vs. viscoelastic-aneurysm: 1D comparison
% Physiological smooth-coefficient model (Sec. 3.2/App. B): Eq.(17) is
% solved over the ENTIRE domain with smooth coefficients, no internal
% jump conditions. Three cases share the identical domain length, inlet
% pulse, boundary conditions, numerical scheme, mesh, and time step:
%   Case 1: fully healthy elastic artery (no aneurysm)
%   Case 2: artery with an ELASTIC aneurysm (wall viscosity = 0)
%   Case 3: artery with a VISCOELASTIC (Kelvin-Voigt) aneurysm
% Cases 2 and 3 share IDENTICAL aneurysm geometry (r_a,h_a,E_a); they
% differ ONLY in wall viscosity eta_a, isolating its effect.

clear; clc; close all; 

%% SI parameters and characteristic scales (Table 3)
rho = 1056; mu = 3.5e-3;
L = 0.150; T = 0.050;
A_star = 3.06e-7;
Q_star = L*A_star/T;

%% Healthy endpoint values (shared by all three cases)
r_h = 1.425e-3; h_h = 2.5e-4;
k1 = 2.00e6; k2 = -2.253e3; k3 = 8.65e4;      % Olufsen empirical relation
E_h = (r_h/h_h)*(k1*exp(k2*r_h)+k3);
Aref_h = pi*r_h^2;
alpha_h = E_h*h_h/(2*pi*r_h^3);
c_h = sqrt(alpha_h*Aref_h/rho); C_h = c_h*T/L;

%% Table 3 aneurysm cross-section/material (defines cases 2 and 3)
r_a = 2.50e-3; h_a = 1.50e-4; E_a = 1.00e6; eta_a = 2.00e3;

%% Shared geometry: aneurysm location, neck transition width, mesh
z_a1 = 0.0725; z_a2 = 0.0775;                 % m (Table 3)
x1 = z_a1/L; x2 = z_a2/L;
w_neck = 0.25*(x2-x1);   % modelling choice (paper leaves the width unspecified)
N_requested = 500;
N_min = ceil(10/(2*w_neck));                  % >=10 grid points across the transition
N = max(N_requested,N_min);
dx = 1/N; Nu = N+1;
x_plot = (0:N)'*dx;

%% Shared time grid, CFL-adaptive (Eq. 39). C_h bounds the wave speed in
% ALL three cases (verified: primitive-based interpolation keeps c(z)
% monotonic between c_h and c_a, with c_h the larger of the two), so one
% shared dtau/dx is valid and CFL-safe for all three cases identically.
t_final = 2.0; dt_phys_requested = 2.0e-5;
tau_final = t_final/T;
cfl_safety = 0.9;
dtau_cfl_max = dx/C_h;
dtau = min(dt_phys_requested/T, cfl_safety*dtau_cfl_max);
Nt = round(tau_final/dtau); dtau = tau_final/Nt;
tau = (0:Nt)*dtau; t = T*tau;
fprintf('Shared mesh/time step (identical for all three cases):\n');
fprintf('N=%d, dx=%.4e, dtau=%.4e, Nt=%d, CFL=%.4f\n',N,dx,dtau,Nt,C_h*dtau/dx);

%% Inlet forcing: triple-sech pulse (Table 1), C-infinity activation (Eq.11-13)
P_P=25*133.322; P_T=15*133.322; P_D=12.5*133.322;
t_P=0.05; t_T=0.20; t_D=0.38; w_P=0.061; w_T=0.076; w_D=0.091;
t_act=0.05; t_delay=0.20;
p_raw = @(tt) P_P./cosh((tt-t_delay-t_P)/w_P) + P_T./cosh((tt-t_delay-t_T)/w_T) + ...
              P_D./cosh((tt-t_delay-t_D)/w_D);
S_inf = @(tt) activation(tt,t_act);
p_in  = @(tt) S_inf(tt).*p_raw(tt);
f_2   = @(ta) C_h*p_in(T*ta)/(alpha_h*A_star);      % = q(0,tau), dimensionless inlet BC
Qtilde_in = @(tt) (c_h/alpha_h)*p_in(tt);           % Sec. 3.4: Q_tilde(zin,t)~(c_h/alpha_h)p_in(t)

%% ================= INLET FUNCTION CHECK =================
[p_in_max,idx_pmax] = max(p_in(t));
[Qin_max,idx_Qmax]  = max(abs(Qtilde_in(t)));
[f2_max,idx_f2max]  = max(abs(f_2(tau)));

fprintf('\n=== Inlet function check ===\n');
fprintf('max p_in        = %.6e Pa    at t=%.6e s (tau=%.4f)\n',p_in_max,t(idx_pmax),tau(idx_pmax));
fprintf('max |Qtilde_in| = %.6e m^3/s at t=%.6e s (tau=%.4f)\n',Qin_max,t(idx_Qmax),tau(idx_Qmax));
fprintf('max |f_2|       = %.6f (dimensionless) at tau=%.4f (t=%.6e s)\n',f2_max,tau(idx_f2max),t(idx_f2max));

f2_0 = f_2(0);
fprintf('\nf_2(0) = %.3e (should be ~0)\n',f2_0);
fprintf('f_2''(0), forward-difference estimate at shrinking eps (should -> 0):\n');
for eps_test = [10*dtau, dtau, dtau/10, dtau/100]
    df2_est = (f_2(eps_test)-f_2(0))/eps_test;
    fprintf('  eps=%.3e:  f_2''(0) ~ %.3e\n',eps_test,df2_est);
end

figure('Color','w'); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(t,p_in(t),'LineWidth',1.5); grid on; box on;
xlabel('$t\;[\mathrm{s}]$','Interpreter','latex'); ylabel('$p_{in}(T\tau)\;[\mathrm{Pa}]$','Interpreter','latex');
title('Inlet pressure $p_{in}(T\tau)$','Interpreter','latex');
nexttile; plot(t,Qtilde_in(t),'LineWidth',1.5); grid on; box on;
xlabel('$t\;[\mathrm{s}]$','Interpreter','latex'); ylabel('$\tilde{Q}_{in}(T\tau)\;[\mathrm{m^3/s}]$','Interpreter','latex');
title('Inlet flow-rate perturbation $\tilde{Q}_{in}(T\tau)$','Interpreter','latex');
nexttile; plot(tau,f_2(tau),'LineWidth',1.5); grid on; box on;
xlabel('$\tau$','Interpreter','latex'); ylabel('$f_2(\tau)$','Interpreter','latex');
title('Dimensionless inlet BC $f_2(\tau)=q(0,\tau)$','Interpreter','latex');
nexttile; plot(t,S_inf(t),'LineWidth',1.5); grid on; box on; ylim([-0.05 1.05]);
xlabel('$t\;[\mathrm{s}]$','Interpreter','latex'); ylabel('$S_\infty(t)$','Interpreter','latex');
title('Activation function $S_\infty(t)$','Interpreter','latex');
sgtitle('Inlet function check');

%% ============ Activation-time sensitivity: inlet waveform ======= %%

t_act_values = [0.05 0.10 0.20 0.24 0.50];
figure('Color','w'); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% Activation functions
nexttile; hold on;
for ia = 1:numel(t_act_values)
    plot(t,activation(t,t_act_values(ia)),'LineWidth',1.3, 'DisplayName',sprintf('$t_{act}=%.2f$ s',t_act_values(ia)));
end
grid on; box on; xlabel('$t$ [s]','Interpreter','latex'); ylabel('$S_\infty(t)$','Interpreter','latex'); title('Activation function','Interpreter','latex');
ylim([-0.05 1.05]); legend('Interpreter','latex','Location','best');

% Inlet pressure
nexttile; hold on;
for ia = 1:numel(t_act_values)
    tact_test = t_act_values(ia); p_test = activation(t,tact_test).*p_raw(t);
    plot(t,p_test/133.322,'LineWidth',1.3,'DisplayName',sprintf('$t_{act}=%.2f$ s',tact_test));
end
grid on; box on; xlabel('$t$ [s]','Interpreter','latex');ylabel('$p_{in}$ [mmHg]','Interpreter','latex'); title('Inlet pressure','Interpreter','latex');
legend('Interpreter','latex','Location','best');sgtitle('Sensitivity to inlet activation time');

%% =========== Activation-time sensitivity: frequency content ======== %%
t_fft_max = 1.0; fft_mask = t <= t_fft_max; t_fft = t(fft_mask); dt_fft = t_fft(2)-t_fft(1); n_fft = numel(t_fft);
figure('Color','w'); hold on;
for ia = 1:numel(t_act_values)
    tact_test = t_act_values(ia); p_test = activation(t_fft,tact_test).*p_raw(t_fft);p_test = p_test - mean(p_test);
    Y = fft(p_test); P2 = abs(Y/n_fft);P1 = P2(1:floor(n_fft/2)+1); P1(2:end-1) = 2*P1(2:end-1);f = (0:floor(n_fft/2))/(n_fft*dt_fft);
    % Avoid plotting zero values on logarithmic scale
    P1(P1 <= 0) = NaN; semilogy(f,P1,'LineWidth',1.3,'DisplayName',sprintf('$t_{act}=%.2f$ s',tact_test));
end
xlim([0 15]); grid on; box on; xlabel('frequency [Hz]'); ylabel('amplitude');

title('Inlet-pressure frequency content'); legend('Interpreter','latex','Location','best');

%% ================= THREE-CASE COMPARISON =================
cases = struct( ...
    'name',  {'Case 1: Fully healthy','Case 2: Elastic aneurysm','Case 3: Viscoelastic aneurysm'}, ...
    'r_a',   {r_h, r_a, r_a}, ...
    'h_a',   {h_h, h_a, h_a}, ...
    'E_a',   {E_h, E_a, E_a}, ...
    'eta_a', {0,   0,   eta_a});

results(3) = struct('q',[],'a',[],'p',[],'Q',[]);
for ic = 1:3
    cs = cases(ic);
    [q_c,a_c,p_c,Q_c,C_case] = simulate_case(cs.r_a,cs.h_a,cs.E_a,cs.eta_a, ...
        rho,mu,L,T,A_star,Q_star, r_h,h_h,E_h,C_h,alpha_h, ...
        x1,x2,w_neck, N,dx,Nu,x_plot, dtau,Nt,tau, f_2);
    results(ic).q = q_c; results(ic).a = a_c; results(ic).p = p_c; results(ic).Q = Q_c;
    fprintf('\n%s: CFL=%.4f, peak q_out/q_in=%.6f\n', ...
            cs.name, max(C_case)*dtau/dx, max(q_c(end,:))/max(q_c(1,:)));
end

%% ============== Plots: probe waveforms and space-time maps, per case ====%%
x_probe = [0.25, 0.50, 0.75];
probe_names = {'healthy proximal','aneurysm region','healthy distal'};
tau_plot_max = 22;    % restrict to where the pulse is actually active
tmask = tau <= tau_plot_max;
for ic = 1:3
    cs = cases(ic);
    q_c = results(ic).q; a_c = results(ic).a; p_c = results(ic).p;

    figure('Color','w'); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile; hold on;
    for kp = 1:numel(x_probe)
        [~,jp] = min(abs(x_plot-x_probe(kp)));
        plot(tau,q_c(jp,:),'LineWidth',1.3);
    end
    xlim([0 tau_plot_max]); grid on; box on; xlabel('\tau'); ylabel('q');
    legend(probe_names,'Location','best'); title('q(\tau)');
    nexttile; hold on;
    for kp = 1:numel(x_probe)
        [~,jp] = min(abs(x_plot-x_probe(kp)));
        plot(tau,a_c(jp,:),'LineWidth',1.3);
    end
    xlim([0 tau_plot_max]); grid on; box on; xlabel('\tau'); ylabel('a');
    legend(probe_names,'Location','best'); title('a(\tau)');
    nexttile; hold on;
    for kp = 1:numel(x_probe)
        [~,jp] = min(abs(x_plot-x_probe(kp)));
        plot(tau,p_c(jp,:),'LineWidth',1.3);
    end
    xlim([0 tau_plot_max]); grid on; box on; xlabel('\tau'); ylabel('p  [Pa]');
    legend(probe_names,'Location','best'); title('p(\tau)');
    sgtitle(sprintf('%s: probe waveforms',cs.name));

    figure('Color','w'); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile;
    imagesc(x_plot,tau(tmask),q_c(:,tmask)'); axis xy; colormap(gca,'jet'); colorbar;
    % clim(max(abs(q_c(:)))*[-1 1]); hold on;
    xline(x1,'k--','LineWidth',1.3); xline(x2,'k--','LineWidth',1.3);
    xlabel('$x$','Interpreter','latex'); ylabel('$\tau$','Interpreter','latex'); title('$q(x,\tau)$','Interpreter','latex');
    nexttile;
    imagesc(x_plot,tau(tmask),a_c(:,tmask)'); axis xy; colormap(gca,'jet'); colorbar;
    % clim(max(abs(a_c(:)))*[-1 1]); hold on;
    xline(x1,'k--','LineWidth',1.3); xline(x2,'k--','LineWidth',1.3);
    xlabel('$x$','Interpreter','latex'); ylabel('$\tau$','Interpreter','latex'); title('$a(x,\tau)$','Interpreter','latex');
    nexttile;
    imagesc(x_plot,tau(tmask),p_c(:,tmask)'); axis xy; colormap(gca,'jet'); colorbar;
    % clim(max(abs(p_c(:)))*[-1 1]); hold on;
    xline(x1,'k--','LineWidth',1.3); xline(x2,'k--','LineWidth',1.3);
    xlabel('$x$','Interpreter','latex'); ylabel('$\tau$','Interpreter','latex'); title('$p(x,\tau)\;[\mathrm{Pa}]$','Interpreter','latex');
    sgtitle(sprintf('%s: space-time maps (aneurysm location marked)',cs.name));
end
%% ========== Spatial snapshots: flow-wave propagation ======%%

snapshot_tau = [0.10 0.20 0.30 0.40];
% snapshot_tau = [0.2 0.4 5.0 5.4];

figure('Color','w'); tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for ic = 1:3
    cs = cases(ic); q_c = results(ic).q;
    nexttile; hold on;
    % Flow snapshots
    for k = 1:numel(snapshot_tau)
        [~,j] = min(abs(tau-snapshot_tau(k))); plot(x_plot,q_c(:,j),'LineWidth',1.4, 'DisplayName',sprintf('t = %.3f s',t(j)));
    end
    % Mark aneurysm region for Cases 2 and 3
    if ic > 1
        yl = ylim; hA = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.90 0.90 0.90],'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
        uistack(hA,'bottom'); xline(x1,'k--','LineWidth',1.0, 'HandleVisibility','off'); xline(x2,'k--','LineWidth',1.0, 'HandleVisibility','off');
    end
    grid on; box on; xlim([0 1]); xlabel('$x$','Interpreter','latex'); ylabel('$q$','Interpreter','latex'); title(cs.name,'Interpreter','latex');
end

% Common legend: only the four time snapshots
lgd = legend('Location','southoutside', 'Orientation','horizontal'); lgd.Layout.Tile = 'south';

sgtitle({'Flow-wave propagation: spatial snapshots', 'Shaded region indicates the aneurysm location in Cases 2 and 3'});

%% ========= Spatial snapshots: zoom around aneurysm ========= %%

% snapshot_tau = [0.20 0.30 0.40 0.50 0.60];
snapshot_tau = [0.10 0.20 0.30 0.35 0.40 0.50];

figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for ic = 2:3
    cs = cases(ic); q_c = results(ic).q; nexttile; hold on;
    for k = 1:numel(snapshot_tau)
        [~,j] = min(abs(tau-snapshot_tau(k))); plot(x_plot,q_c(:,j),'LineWidth',1.4, 'DisplayName',sprintf('t = %.3f s',t(j)));
    end

    yl = ylim; hA = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.90 0.90 0.90], 'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
    uistack(hA,'bottom'); xline(x1,'k--','LineWidth',1.0,'HandleVisibility','off'); xline(x2,'k--','LineWidth',1.0,'HandleVisibility','off');
    grid on; box on; xlim([0.40 0.60]); xlabel('$x$','Interpreter','latex'); ylabel('$q$','Interpreter','latex'); title(cs.name,'Interpreter','latex');
end

lgd = legend('Location','southoutside', 'Orientation','horizontal'); lgd.Layout.Tile = 'south';

sgtitle({'Flow-wave interaction with the aneurysm', 'Zoomed spatial snapshots for Cases 2 and 3'});
%% ===================== Local functions =====================
function [q,a,p,Q,C] = simulate_case(r_a_eff,h_a_eff,E_a_eff,eta_a_eff, ...
        rho,mu,L,T,A_star,Q_star, r_h,h_h,E_h,C_h,alpha_h, ...
        x1,x2,w_neck, N,dx,Nu,x_plot, dtau,Nt,tau, f_2)
    % Runs ONE case with the SHARED mesh (N,dx), time grid (dtau,Nt,tau),
    % and inlet forcing (f_2) supplied by the caller -- guaranteeing all
    % three cases use identical numerics, differing only in the trial
    % aneurysm cross-section/material (r_a_eff,h_a_eff,E_a_eff,eta_a_eff).
    % r_a_eff=r_h,h_a_eff=h_h,E_a_eff=E_h,eta_a_eff=0 reduces to a fully
    % healthy artery (Case 1); cases 2/3 share r_a_eff,h_a_eff,E_a_eff,
    % differing only in eta_a_eff (elastic vs. viscoelastic).
    phi = neck_localization(x_plot,x1,x2,w_neck);
    r_z = r_h + (r_a_eff-r_h)*phi;
    h_z = h_h + (h_a_eff-h_h)*phi;
    E_z = E_h + (E_a_eff-E_h)*phi;
    Aref_z  = pi*r_z.^2;
    alpha_z = E_z.*h_z./(2*pi*r_z.^3);
    c_z = sqrt(alpha_z.*Aref_z/rho);
    gamma_z = 8*pi*mu./(rho*Aref_z);

    Aref_a_eff = pi*r_a_eff^2;
    beta_a_eff = eta_a_eff*h_a_eff/(2*pi*r_a_eff^3);
    D_a_eff = beta_a_eff*Aref_a_eff/rho;
    D_z = D_a_eff*phi;                 % Eq. B.1 literal, D_h=0
    beta_z = D_z.*rho./Aref_z;         % back-derived, consistent with D(z)
    % D_z = D_h + (D_a_eff-D_h)*phi;   % For general case Dh>0
    % D_h = beta_h*Aref_h/rho;         % if one define a healthy-wall beta_h 

    C = c_z*T/L; G = gamma_z*T; V = D_z*T/L^2;

    A = spalloc(Nu,Nu,3*Nu);
    interior = 2:(Nu-1);
    for j = interior
        lam = V(j)/(2*dtau*dx^2);
        A(j,j-1) = -lam; A(j,j) = 1/dtau^2 + G(j)/(2*dtau) + 2*lam; A(j,j+1) = -lam;
    end
    row_in = 1; A(row_in,row_in) = 1;
    row_out = Nu; A(row_out,row_out) = 1/(2*dtau) + C_h/dx; A(row_out,row_out-1) = -C_h/dx;

    q = zeros(Nu,Nt+1);
    q(row_in,:) = f_2(tau);
    for m = 2:Nt
        qm = q(:,m); qold = q(:,m-1); rhs = zeros(Nu,1);
        rhs(row_in) = f_2(tau(m+1));
        for j = interior
            qxx_m   = (qm(j+1)  -2*qm(j)  +qm(j-1)  )/dx^2;
            qxx_old = (qold(j+1)-2*qold(j)+qold(j-1))/dx^2;
            rhs(j) = (2*qm(j)-qold(j))/dtau^2 + G(j)*qold(j)/(2*dtau) + ...
                     C(j)^2*qxx_m - V(j)*qxx_old/(2*dtau);
        end
        rhs(row_out) = qold(row_out)/(2*dtau);
        q(:,m+1) = A\rhs;
    end

    qx = deriv_1d(q,dx);
    a = zeros(Nu,Nt+1);
    for m = 1:Nt, a(:,m+1) = a(:,m) - 0.5*dtau*(qx(:,m+1)+qx(:,m)); end
    p = alpha_z.*(A_star*a) - beta_z.*(A_star/T).*qx;
    Q = Q_star*q;
end

function phi = neck_localization(x,x1,x2,w)
    rise = theta_smoothstep((x-(x1-w))/(2*w));
    fall = theta_smoothstep((x-(x2-w))/(2*w));
    phi = rise.*(1-fall);
end

function y = theta_smoothstep(s)
    y = zeros(size(s)); mid = s>0 & s<1; u = s(mid);
    y(mid) = exp(-1./u)./(exp(-1./u)+exp(-1./(1-u))); y(s>=1) = 1;
end

function qx = deriv_1d(q,dx)
    qx = zeros(size(q));
    qx(2:end-1,:) = (q(3:end,:)-q(1:end-2,:))/(2*dx);
    qx(1,:)   = (-3*q(1,:)  +4*q(2,:)  -q(3,:)  )/(2*dx);
    qx(end,:) = ( 3*q(end,:)-4*q(end-1,:)+q(end-2,:))/(2*dx);
end

function S = activation(t,t_act)
    S = zeros(size(t)); mid = (t>0)&(t<t_act); u = t(mid)/t_act;
    S(mid) = exp(-1./u)./(exp(-1./u)+exp(-1./(1-u))); S(t>=t_act) = 1;
end
