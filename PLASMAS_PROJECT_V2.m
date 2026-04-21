%%Michael Kuang
%%Physics 525 Project 1 
%%University of Wisconsin - Madison
clear; clc; close all
 
%% Constants & Parameters
e    = 1.60217663e-19; %electron charge
m_e  = 9.1093837e-31; %electron mass
eps0 = 8.85419e-12;  %permitivity of free space constant    
k_B  = 1.380649e-23; %boltzmann constant
n_o = 1e10;  %Background number density 
T_e = 100000;  %Electron temperature 
n_1 = 10; %preturbed number densitry 
a = 0.01; %seperation distance of plates in meters
omega_pdep_sq = (n_o * e^2) / (m_e * eps0); %depressed plasma frequency

%% Dispersion Relations
%initialize the vectors used for calculating quantities
k_perp   = linspace(0.02, 10, 500);

%% FIGURE 2 Warm Dispersion Relations UNCOUPLED
A = omega_pdep_sq / 2; %Dimensionless Parameters  
B = (3 * omega_pdep_sq * k_B * T_e) / (2 * m_e); 
disc = (A .* k_perp).^2 + 4 .* B .* k_perp.^3;  
omega_sq = (A .* k_perp + sqrt(disc)) / 2;         
omega_warm = sqrt(omega_sq);
omega_cold = sqrt((omega_pdep_sq/2) * k_perp);

figure;
hold on;
plot(k_perp, omega_cold, "k-", 'LineWidth',1);
plot(k_perp, omega_warm, "r-", 'LineWidth',1);
title("Warm Plasma Dispersion Relations, T = 100000");
xlabel("wavenumber k_{perp} (1/m)");
ylabel("frequecy \omega (rad/s)");
legend('Cold Plasma', 'Warm Plasma');
hold off

%% FIGURE 3: Coupled Warm Dispersion Relations
C = (3 * k_B * T_e) / m_e;
omega_warm_coupled_1 = zeros(size(k_perp)); 
omega_warm_coupled_2 = zeros(size(k_perp));
omega_cold_coupled_1 = @(k_perp) sqrt(((omega_pdep_sq/2) * k_perp)*(1+exp(-k_perp*a))); %calculation of cold limit
omega_cold_coupled_2 = @(k_perp) sqrt(((omega_pdep_sq/2) * k_perp)*(1-exp(-k_perp*a)));
for i = 1:length(k_perp)
    omega_cold_coupled_values_1(i) = omega_cold_coupled_1(k_perp(i));
    omega_cold_coupled_values_2(i) = omega_cold_coupled_2(k_perp(i));
    
end
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
plot(k_perp, omega_warm_coupled_1,'r-', 'LineWidth', 1.8, 'DisplayName', 'Warm Coupled \omega_+ (high)');
plot(k_perp, omega_warm_coupled_2,'b-', 'LineWidth', 1.8, 'DisplayName', 'Warm Coupled \omega_- (low)');
plot(k_perp, omega_cold,'k:', 'LineWidth', 1.0, 'DisplayName', 'Uncoupled Cold Limit');
title(sprintf('Warm vs Cold Coupled Plasma Dispersion Relations, T_e = %.0f K, a = %.3f m', T_e, a));
xlabel('k_{\perp} (1/m)');
ylabel('\omega (rad/s)');
legend('show','Location', 'northwest');
grid on;
hold off;

%% Figure 4: COUPLED WARM GROUP VELOCITY VS SEPERATION DISTANCE
a_range = linspace(0.01, 10, 400);
k_vals  = [0.5, 1.0];
colors  = {'b', 'r'};

figure;
hold on;
for ki = 1:length(k_vals)
    k = k_vals(ki);
    vg_plus = zeros(size(a_range));
    vg_minus = zeros(size(a_range));
    for j = 1:length(a_range)
        a_j = a_range(j);
        ep  = exp(-k * a_j);
        fac_p = (1 + ep) * omega_pdep_sq * k;
        r_p = roots([2, -fac_p, -fac_p * C * k^2]);
        r_p = r_p(r_p > 0);
        w_p = sqrt(max(r_p));
        dFp_dw = -(2*C*k^2 / w_p^3)*(1+ep) - 4*w_p/(omega_pdep_sq*k);
        dFp_dk =  (2*C*k   / w_p^2)*(1+ep) - a_j*ep*(1 + C*k^2/w_p^2) + 2*w_p^2/(omega_pdep_sq*k^2);
        vg_plus(j) = -dFp_dk / dFp_dw;
        fac_m = (1 - ep) * omega_pdep_sq * k;
        if fac_m > 0
            r_m = roots([2, -fac_m, -fac_m * C * k^2]);
            r_m = r_m(r_m > 0);
            if ~isempty(r_m)
                w_m = sqrt(max(r_m));
                dFm_dw = -(2*C*k^2 / w_m^3)*(1-ep) - 4*w_m/(omega_pdep_sq*k);
                dFm_dk = (2*C*k   / w_m^2)*(1-ep) + a_j*ep*(1 + C*k^2/w_m^2) + 2*w_m^2/(omega_pdep_sq*k^2);
                vg_minus(j) = -dFm_dk / dFm_dw;
            else
                vg_minus(j) = NaN;
            end
        else
            vg_minus(j) = NaN;
        end
    end
    plot(a_range, vg_plus,  '-',  'Color', colors{ki}, 'LineWidth', 1.5,'DisplayName', sprintf('\\omega_+ k_\\perp = %.1f m^{-1}', k));
    plot(a_range, vg_minus, '--', 'Color', colors{ki}, 'LineWidth', 1.5,'DisplayName', sprintf('\\omega_-   k_\\perp = %.1f m^{-1}', k));
    B_k = (3 * omega_pdep_sq * k_B * T_e) / (2 * m_e);
    disc_k = (A * k)^2 + 4 * B_k * k^3;
    omega_k = sqrt((A * k + sqrt(disc_k)) / 2);
    dosq_dk = (A + (A^2 * k + 6*B_k*k^2) / sqrt(disc_k)) / 2;
    vg_warm_k = dosq_dk / (2 * omega_k);
    yline(vg_warm_k, ':', 'Color', colors{ki}, 'LineWidth', 1.2,'DisplayName', sprintf('Warm Uncoupled Limit  k_\\perp = %.1f m^{-1}', k));
end
title(sprintf('Warm Plasma Group Velocity vs Separation Distance, T_e = %.0f K', T_e));
xlabel('Separation distance a (m)');
ylabel('v_g (m/s)');
legend('show', 'Location', 'best');
grid on;
hold off;

%% Figure 5: Group velocity vs temperature parameterized by wavenumber
%uncoupled 
T_range = linspace(100, 1e5, 300);
k_array = [0.5, 1.0, 2.0, 3.0, 5.0];
colors_T = {'b-','r-','k-','g-','m-'};
figure;
hold on;
for ki = 1:length(k_array)
    k_fixed = k_array(ki);
    vg_vs_T = zeros(size(T_range));
    for j = 1:length(T_range)
        B_j = (3 * omega_pdep_sq * k_B * T_range(j)) / (2 * m_e);
        disc_j = (A * k_fixed)^2 + 4 * B_j * k_fixed^3;
        omega_j = sqrt((A * k_fixed + sqrt(disc_j)) / 2);
        d_osq_j = (A + (A^2 * k_fixed + 6*B_j*k_fixed^2) / sqrt(disc_j)) / 2;
        vg_vs_T(j) = d_osq_j / (2 * omega_j);
    end
    plot(T_range, vg_vs_T, colors_T{ki}, 'LineWidth', 1.5,'DisplayName', sprintf('k_\\perp = %.1f m^{-1}', k_fixed));
end
title('Warm Plasma Group Velocity vs Temperature (multiple k_{\perp})');
xlabel('T_e (K)');
ylabel('v_g (m/s)');
legend('show');
grid on;
hold off;

%% Figure 6: Warm plasma uncoupled group velocity vs wavenumber 
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