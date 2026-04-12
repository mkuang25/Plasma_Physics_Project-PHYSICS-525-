%Michael Kuang
%2D Plasma Oscillations
%Dispersion relations
% -------------------------------------------------------------------------
clear; clc; close all;

%%Constants & Parameters
e    = 1.60217663e-19; %electron charge
m_e  = 9.1093837e-31; %electron mass
eps0 = 8.85419e-12;  %permitivity of free space constant    
k_B  = 1.380649e-23; %boltzmann constant
n_o = 1e10;  %Background number density 
T_e = 100000;  %Electron temperature 
n_1 = 10; %preturbed number densitry 
a = 0.01; %seperation distance of plates in meters
omega_pdep_sq = (n_o * e^2) / (m_e * eps0);


%% Dispersion Relation
%Cold Dispersion Relations 
k_perp   = linspace(0.02, 10, 500);
omega_cold = @(k_perp) sqrt((omega_pdep_sq/2) * k_perp);
omega_cold_coupled_1 = @(k_perp) sqrt(((omega_pdep_sq/2) * k_perp)*(1+exp(-k_perp*a)));
omega_cold_coupled_2 = @(k_perp) sqrt(((omega_pdep_sq/2) * k_perp)*(1-exp(-k_perp*a)));
for i = 1:length(k_perp)
    omega_cold_values(i) = omega_cold(k_perp(i));
    omega_cold_coupled_values_1(i) = omega_cold_coupled_1(k_perp(i));
    omega_cold_coupled_values_2(i) = omega_cold_coupled_2(k_perp(i));
    
end

%Warm Dispersion Relations
A = omega_pdep_sq / 2;   
B = (3 * omega_pdep_sq * k_B * T_e) / (2 * m_e); 
disc = (A .* k_perp).^2 + 4 .* B .* k_perp.^3;  
omega_sq = (A .* k_perp + sqrt(disc)) / 2;         
omega_warm = sqrt(omega_sq);

%% Group Velocity
% --- Figure 3: Group velocity vs wavenumber (warm plasma) ---
d_omega_sq_dk = (A + (A^2 .* k_perp + 6.*B.*k_perp.^2) ./ sqrt(disc)) / 2;
vg_warm = d_omega_sq_dk ./ (2 .* omega_warm);
vg_cold = 0.5 *(sqrt(omega_pdep_sq/2)) .* (1./(sqrt(k_perp)));

figure;
hold on
plot(k_perp, vg_warm, 'r-', 'LineWidth', 1.5);
plot(k_perp, vg_cold, 'b-', 'LineWidth', 1.5);
title('Warm Plasma Group Velocity vs Wavenumber, T = 100000');
xlabel('k_{\perp} (1/m)');
ylabel('v_g (m/s)');
legend('Warm Plasma', 'Cold Plasma');
grid on;
hold off;

% --- Figure 4: Group velocity vs Temperature (multiple k_perp) ---
T_range   = linspace(100, 1e5, 300);
k_array   = [0.5, 1.0, 2.0, 3.0, 5.0];
colors_T  = {'b-','r-','k-','g-','m-'};

figure;
hold on;
for ki = 1:length(k_array)
    k_fixed  = k_array(ki);
    vg_vs_T  = zeros(size(T_range));
    for j = 1:length(T_range)
        B_j        = (3 * omega_pdep_sq * k_B * T_range(j)) / (2 * m_e);
        disc_j     = (A * k_fixed)^2 + 4 * B_j * k_fixed^3;
        omega_j    = sqrt((A * k_fixed + sqrt(disc_j)) / 2);
        d_osq_j    = (A + (A^2 * k_fixed + 6*B_j*k_fixed^2) / sqrt(disc_j)) / 2;
        vg_vs_T(j) = d_osq_j / (2 * omega_j);
    end
    plot(T_range, vg_vs_T, colors_T{ki}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('k_\\perp = %.1f m^{-1}', k_fixed));
end
title('Warm Plasma Group Velocity vs Temperature (multiple k_{\perp})');
xlabel('T_e (K)');
ylabel('v_g (m/s)');
legend('show');
grid on;
hold off;

% --- Figure 5: Group velocity vs separation distance a (multiple k_perp) ---
a_range   = linspace(0.01, 1, 400);
k_array2  = [0.5, 1.0];
colors_p  = {'b-','r-'};

figure;
hold on;
for ki = 1:length(k_array2)
    k_fixed2 = k_array2(ki);
    vg_mode1 = zeros(size(a_range));
    vg_mode2 = zeros(size(a_range));
    for j = 1:length(a_range)
        a_j          = a_range(j);
        ep           = exp(-k_fixed2 * a_j);
        omega1       = sqrt(A * k_fixed2 * (1 + ep));
        dg1_dk       = A * (1 + ep*(1 - k_fixed2*a_j));
        vg_mode1(j)  = dg1_dk / (2 * omega1);
        omega2       = sqrt(A * k_fixed2 * (1 - ep));
        dg2_dk       = A * (1 + ep*(k_fixed2*a_j - 1));
        vg_mode2(j)  = dg2_dk / (2 * omega2);
    end
    plot(a_range, vg_mode1, colors_p{ki}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\omega_+  k_\\perp = %.1f', k_fixed2));
    plot(a_range, vg_mode2, colors_p{ki}, 'LineWidth', 1.0, 'LineStyle', '--', ...
        'DisplayName', sprintf('\\omega_-  k_\\perp = %.1f', k_fixed2));
end
title('Coupled Cold Plasma Group Velocities vs Separation (multiple k_{\perp})');
xlabel('Separation distance a (m)');
ylabel('v_g (m/s)');
legend('show', 'Location', 'best');
grid on;
hold off;

% -- Figure 6: Warm Coupled Plasmas 
C = (3 * k_B * T_e) / m_e;
omega_warm_coupled_1 = zeros(size(k_perp));
omega_warm_coupled_2 = zeros(size(k_perp));
for i = 1:length(k_perp)
    k = k_perp(i);
    ep = exp(-k * a);

    fac1 = (1 + ep) * omega_pdep_sq * k;
    coeffs1 = [2, -fac1, -fac1 * C * k^2];
    r1 = roots(coeffs1);
    r1 = r1(r1 > 0);
    omega_warm_coupled_1(i) = sqrt(max(r1));

    fac2 = (1 - ep) * omega_pdep_sq * k;
    coeffs2 = [2, -fac2, -fac2 * C * k^2];
    r2 = roots(coeffs2);
    r2 = r2(r2 > 0);
    if ~isempty(r2)
        omega_warm_coupled_2(i) = sqrt(max(r2));
    else
        omega_warm_coupled_2(i) = NaN;
    end
end

figure;
hold on;
plot(k_perp, omega_cold_coupled_values_1, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Cold Coupled \omega_+');
plot(k_perp, omega_cold_coupled_values_2, 'b--', 'LineWidth', 1.2, 'DisplayName', 'Cold Coupled \omega_-');
plot(k_perp, omega_warm_coupled_1,        'r-',  'LineWidth', 1.8, 'DisplayName', 'Warm Coupled \omega_+ (high)');
plot(k_perp, omega_warm_coupled_2,        'b-',  'LineWidth', 1.8, 'DisplayName', 'Warm Coupled \omega_- (low)');
plot(k_perp, omega_cold_values,           'k:',  'LineWidth', 1.0, 'DisplayName', 'Uncoupled Cold Limit');
title(sprintf('Warm vs Cold Coupled Plasma Dispersion Relations, T_e = %.0f K, a = %.3f m', T_e, a));
xlabel('k_{\perp} (1/m)');
ylabel('\omega (rad/s)');
legend('show', 'Location', 'northwest');
grid on;
hold off;

%% Figure 7: Warm Coupled Group Velocities vs Separation Distance
C        = (3 * k_B * T_e) / m_e;
k_array3 = [0.5, 1.0];
colors_w = {'b-','r-'};

figure;
hold on;
for ki = 1:length(k_array3)
    k     = k_array3(ki);
    vg_w1 = zeros(size(a_range));
    vg_w2 = zeros(size(a_range));

    for j = 1:length(a_range)
        a_j = a_range(j);
        ep  = exp(-k * a_j);

        % High-frequency mode (+): fac1 = omega_p^2 * k * (1 + exp(-k*a))
        fac1     = (1 + ep) * omega_pdep_sq * k;
        dfac1_dk = omega_pdep_sq * (1 + ep*(1 - k*a_j));       % d(fac1)/dk
        r1       = roots([2, -fac1, -fac1*C*k^2]);
        r1       = r1(r1 > 0);
        x1       = max(r1);
        omega1   = sqrt(x1);
        dx1_dk   = (dfac1_dk*(x1 + C*k^2) + 2*C*k*fac1) / (4*x1 - fac1);
        vg_w1(j) = dx1_dk / (2 * omega1);

        % Low-frequency mode (-): fac2 = omega_p^2 * k * (1 - exp(-k*a))
        fac2     = (1 - ep) * omega_pdep_sq * k;
        dfac2_dk = omega_pdep_sq * (1 + ep*(k*a_j - 1));       % d(fac2)/dk
        r2       = roots([2, -fac2, -fac2*C*k^2]);
        r2       = r2(r2 > 0);
        if ~isempty(r2)
            x2       = max(r2);
            omega2   = sqrt(x2);
            dx2_dk   = (dfac2_dk*(x2 + C*k^2) + 2*C*k*fac2) / (4*x2 - fac2);
            vg_w2(j) = dx2_dk / (2 * omega2);
        else
            vg_w2(j) = NaN;
        end
    end

    plot(a_range, vg_w1, colors_w{ki}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\omega_+  k_\\perp = %.1f', k));
    plot(a_range, vg_w2, colors_w{ki}, 'LineWidth', 1.0, 'LineStyle', '--', ...
        'DisplayName', sprintf('\\omega_-  k_\\perp = %.1f', k));
end
title(sprintf('Warm Coupled Plasma Group Velocities vs Separation, T_e = %.0f K', T_e));
xlabel('Separation distance a (m)');
ylabel('v_g (m/s)');
legend('show', 'Location', 'best');
grid on;
hold off;

%Fig 1: Dispersion Relations for Cold Plasmas in Coupled modes
figure;
hold on;, 
plot(k_perp, omega_cold_values, "k-", 'LineWidth',1);
plot(k_perp, omega_cold_coupled_values_1, "r-", 'LineWidth',1);
plot(k_perp, omega_cold_coupled_values_2, "b-"','LineWidth',1);
title("Coupled Cold Plasma Dispersion Relations a = 0.01");
xlabel("wavenumber k_{perp} (1/m)");
ylabel("frequecy \omega (rad/s)");
legend('Uncoupled Cold Plasma Limit', 'Coupled Cold Plasma Freq 1', 'Coupled Cold Plasma Freq 2');
hold off

%Fig 2: Dispersion Relations in Warm Plasma vs Cold Plasma UNCOUPLED
figure;
hold on;
plot(k_perp, omega_cold_values, "k-", 'LineWidth',1);
plot(k_perp, omega_warm, "r-", 'LineWidth',1);
title("Temperature Dependent Plasma Dispersion Relations, T = 100000");
xlabel("wavenumber k_{perp} (1/m)");
ylabel("frequecy \omega (rad/s)");
legend('Cold Plasma', 'Warm Plasma');
hold off