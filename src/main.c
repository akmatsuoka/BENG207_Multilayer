/* main.c
 *
 * BENG207_1: Multilayer Acoustic Transfer-Matrix Model (Extended)
 * with Kelvin-Voigt viscoelastic couplant
 *
 * Stack (bottom to top):
 *   1. LiNbO3 transducer (semi-infinite source)
 *   2. Ultrasound gel couplant (Kelvin-Voigt: E* = E + iωη)
 *   3. Borosilicate glass superstrate (e.g., 1 mm)
 *   4. LiNbO3 coverslip (e.g., 150 um)
 *   5. Water fluid channel (terminating half-space)
 *
 * See README file for detail.
 *
 * CROSS-PLATFORM: compiles on Linux, Mac, and Windows (MSVC & MinGW)
 * from the same source file. No <complex.h>, no POSIX dependencies.
 *
 * Compile (Linux/Mac):
 *   gcc -std=c99 -O2 -o tm_model main.c -lm
 *
 * Compile (Windows MSVC — Developer Command Prompt is needed):
 *   cl /O2 main.c /link /out:tm_model.exe
 *
 * Compile (Windows MinGW):
 *   gcc -O2 -o tm_model.exe main.c -lm
 *
 * Run:
 *   ./tm_model          (Linux/Mac)
 *   tm_model.exe        (Windows)
 *   CSVs appear in PLOT_CSV/ subfolder.
 * 
 * Github sync instruction
 * 1. Go to your local directory (e.g., /Projects/BENG207_1/)
 * 2. git int
 * 3. git add .
 * 4. git commit -m "Comment.."
 * 5. git branch -M main
 * 6. git push -u origin main
 * 
 *  ====Other useful tips===================================================
 * Please do your work either on your Windows laptop or on our Leichtag Linus 
 * server. 
 * Please do not change pointer codes.
 *  =========================================================================
 */

 * 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================================================
 * Platform-specific directory creation to avoid the "WINDOWS" issues
 * ============================================================ */
#include <errno.h>
#ifdef _WIN32
  #include <direct.h>
  #define MKDIR(path) _mkdir(path)
#else
  #include <sys/stat.h>
  #define MKDIR(path) mkdir(path, 0755)
#endif

/* ============================================================
 * Complex number helpers (no <complex.h> needed — MSVC lacks it)
 *
 * We represent complex numbers as a struct with .re and .im
 * Without these, we would be suffering a lot!!
 * ============================================================ */

typedef struct { double re; double im; } Cpx;

static Cpx cpx(double re, double im)
{
    Cpx z;
    z.re = re;
    z.im = im;
    return z;
}

static Cpx cpx_add(Cpx a, Cpx b)
{
    return cpx(a.re + b.re, a.im + b.im);
}

static Cpx cpx_sub(Cpx a, Cpx b)
{
    return cpx(a.re - b.re, a.im - b.im);
}

static Cpx cpx_mul(Cpx a, Cpx b)
{
    return cpx(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}

static Cpx cpx_div(Cpx a, Cpx b)
{
    double d = b.re*b.re + b.im*b.im;
    return cpx((a.re*b.re + a.im*b.im) / d,
               (a.im*b.re - a.re*b.im) / d);
}

static Cpx cpx_scale(double s, Cpx a)
{
    return cpx(s * a.re, s * a.im);
}

static double cpx_abs(Cpx a)
{
    return sqrt(a.re*a.re + a.im*a.im);
}

static double cpx_arg(Cpx a)
{
    return atan2(a.im, a.re);
}

static Cpx cpx_exp(Cpx a)
{
    double e = exp(a.re);
    return cpx(e * cos(a.im), e * sin(a.im));
}

static Cpx cpx_sqrt(Cpx a)
{
    double r = sqrt(cpx_abs(a));
    double theta = cpx_arg(a) / 2.0;
    return cpx(r * cos(theta), r * sin(theta));
}

/* sin(z) = (exp(iz) - exp(-iz)) / (2i) */
static Cpx cpx_sin(Cpx a)
{
    return cpx(sin(a.re) * cosh(a.im), cos(a.re) * sinh(a.im));
}

/* cos(z) = (exp(iz) + exp(-iz)) / 2 */
static Cpx cpx_cos(Cpx a)
{
    return cpx(cos(a.re) * cosh(a.im), -sin(a.re) * sinh(a.im));
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Output directory — relative to executable */
#define OUTPUT_DIR "PLOT_CSV"

/* ============================================================
 * Material parameters
 * ============================================================ */

static const double rho_1     = 4647.0;    /* LiNbO3 kg/m^3 */
static const double c_1_long  = 6570.0;    /* LiNbO3 m/s */

static const double rho_2   = 1020.0;      /* couplant kg/m^3 */
static const double c_2     = 1500.0;      /* couplant m/s */
static const double E_2_0   = 5.0e3;       /* Pa (fresh gel) */
static const double eta_2_0 = 2.0;         /* Pa-s (fresh gel) */

static const double rho_3 = 2230.0;        /* glass kg/m^3 */
static const double c_3   = 5640.0;        /* glass m/s */
static const double h_3   = 1.0e-3;        /* 1 mm */

static const double rho_4 = 4647.0;        /* coverslip kg/m^3 */
static const double c_4   = 6570.0;        /* coverslip m/s */
static const double h_4   = 150.0e-6;      /* 150 um */

static const double rho_w = 1000.0;        /* water kg/m^3 */
static const double c_w   = 1483.0;        /* water m/s */

static const double f_center = 20.0e6;     /* 20 MHz */

static double Z_1, Z_3, Z_4, Z_w;

/* ============================================================
 * 2x2 complex matrix
 * ============================================================ */

typedef struct { Cpx m[2][2]; } Mat2x2;

static Mat2x2 mat_mul(Mat2x2 A, Mat2x2 B)
{
    Mat2x2 C;
    int i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            C.m[i][j] = cpx_add(cpx_mul(A.m[i][0], B.m[0][j]),
                                cpx_mul(A.m[i][1], B.m[1][j]));
    return C;
}

/* ============================================================
 * Kelvin-Voigt complex wave speed
 * c*(omega) = sqrt( (E + i*omega*eta) / rho )
 * ============================================================ */

static Cpx kv_wave_speed(double omega, double E, double eta, double rho)
{
    Cpx E_star = cpx(E, omega * eta);
    return cpx_sqrt(cpx_scale(1.0 / rho, E_star));
}

/* ============================================================
 * Transfer matrix for a single layer
 * ============================================================ */

static Mat2x2 layer_matrix(double omega, double rho, Cpx c_cmplx, double h)
{
    Cpx k  = cpx_div(cpx(omega, 0), c_cmplx);
    Cpx Z  = cpx_scale(rho, c_cmplx);
    Cpx kh = cpx_scale(h, k);
    Cpx cos_kh = cpx_cos(kh);
    Cpx sin_kh = cpx_sin(kh);
    Cpx I_unit = cpx(0, 1);
    Mat2x2 M;

    M.m[0][0] = cos_kh;
    M.m[0][1] = cpx_div(cpx_mul(I_unit, sin_kh), Z);
    M.m[1][0] = cpx_mul(cpx_mul(I_unit, Z), sin_kh);
    M.m[1][1] = cos_kh;
    return M;
}

/* ============================================================
 * Transmission coefficient — generalized (variable h3, h4)
 * ============================================================ */

static Cpx compute_T_full(double freq, double h2, double E2,
                           double eta2, double h3_var, double h4_var)
{
    double omega = 2.0 * M_PI * freq;
    Cpx c2;
    Mat2x2 M2, M3, M4, Mt;
    Cpx denom;

    c2 = kv_wave_speed(omega, E2, eta2, rho_2);
    M2 = layer_matrix(omega, rho_2, c2, h2);
    M3 = layer_matrix(omega, rho_3, cpx(c_3, 0), h3_var);
    M4 = layer_matrix(omega, rho_4, cpx(c_4, 0), h4_var);

    Mt = mat_mul(M2, M3);
    Mt = mat_mul(Mt, M4);

    /* T = 2 / (M11 + M12*Zw + M21/Z1 + M22*Zw/Z1) */
    denom = cpx_add(
                cpx_add(Mt.m[0][0],
                        cpx_scale(Z_w, Mt.m[0][1])),
                cpx_add(cpx_scale(1.0/Z_1, Mt.m[1][0]),
                        cpx_scale(Z_w/Z_1, Mt.m[1][1])));

    return cpx_div(cpx(2.0, 0), denom);
}

static Cpx compute_T(double freq, double h2, double E2, double eta2)
{
    return compute_T_full(freq, h2, E2, eta2, h_3, h_4);
}

/* ============================================================
 * Utilities
 * ============================================================ */

static void linspace(double start, double end, int n, double *out)
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = start + (end - start) * i / (n - 1);
}

static void logspace(double start_exp, double end_exp, int n, double *out)
{
    int i;
    for (i = 0; i < n; i++) {
        double e = start_exp + (end_exp - start_exp) * i / (n - 1);
        out[i] = pow(10.0, e);
    }
}

static double dmax(double a, double b)
{
    return (a > b) ? a : b;
}

/* Build file path: OUTPUT_DIR/filename */
static void build_path(char *buf, int bufsize, const char *filename)
{
#ifdef _WIN32
    _snprintf(buf, bufsize, "%s\\%s", OUTPUT_DIR, filename);
#else
    snprintf(buf, bufsize, "%s/%s", OUTPUT_DIR, filename);
#endif
}

/* ============================================================
 * OUTPUT 1: Transmission spectrum
 * ============================================================ */

static void output_spectrum(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 2000;
    double freqs[2000];
    double h2_vals[] = {5e-6, 15e-6, 30e-6, 50e-6};
    int Nh = 4;
    double E_deg = 50.0e3, eta_deg = 20.0;
    int i, j;

    build_path(fname, sizeof(fname), "fig1_spectrum.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(15.0e6, 25.0e6, Nf, freqs);

    fprintf(fp, "freq_MHz");
    for (j = 0; j < Nh; j++)
        fprintf(fp, ",T_fresh_h%.0fum,T_deg_h%.0fum",
                h2_vals[j]*1e6, h2_vals[j]*1e6);
    fprintf(fp, "\n");

    for (i = 0; i < Nf; i++) {
        fprintf(fp, "%.6f", freqs[i] / 1e6);
        for (j = 0; j < Nh; j++) {
            Cpx Tf = compute_T(freqs[i], h2_vals[j], E_2_0, eta_2_0);
            Cpx Td = compute_T(freqs[i], h2_vals[j], E_deg, eta_deg);
            fprintf(fp, ",%.8f,%.8f", cpx_abs(Tf), cpx_abs(Td));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 2: Sensitivity maps
 * ============================================================ */

static void output_sensitivity(void)
{
    char fname1[256], fname2[256];
    int Nh = 100, Np = 100;
    double h2_arr[100], eta_arr[100], E_arr[100];
    FILE *fp1, *fp2;
    int i, j;

    build_path(fname1, sizeof(fname1), "fig2a_sensitivity_eta.csv");
    build_path(fname2, sizeof(fname2), "fig2b_sensitivity_E.csv");

    linspace(1e-6, 50e-6, Nh, h2_arr);
    logspace(-1, 2, Np, eta_arr);
    logspace(2, 6, Np, E_arr);

    fp1 = fopen(fname1, "w");
    if (!fp1) { perror(fname1); return; }
    fprintf(fp1, "h2_um,eta2_Pas,T_abs\n");
    for (i = 0; i < Np; i++)
        for (j = 0; j < Nh; j++) {
            Cpx T = compute_T(f_center, h2_arr[j], E_2_0, eta_arr[i]);
            fprintf(fp1, "%.4f,%.6f,%.8f\n", h2_arr[j]*1e6, eta_arr[i], cpx_abs(T));
        }
    fclose(fp1);
    printf("  Written: %s\n", fname1);

    fp2 = fopen(fname2, "w");
    if (!fp2) { perror(fname2); return; }
    fprintf(fp2, "h2_um,E2_Pa,T_abs\n");
    for (i = 0; i < Np; i++)
        for (j = 0; j < Nh; j++) {
            Cpx T = compute_T(f_center, h2_arr[j], E_arr[i], eta_2_0);
            fprintf(fp2, "%.4f,%.2f,%.8f\n", h2_arr[j]*1e6, E_arr[i], cpx_abs(T));
        }
    fclose(fp2);
    printf("  Written: %s\n", fname2);
}

/* ============================================================
 * OUTPUT 3: Phase shift during degradation
 * ============================================================ */

static void output_phase_shift(void)
{
    char fname[256];
    FILE *fp;
    int N = 200;
    double E2_arr[200], eta2_arr[200], h2_arr[200];
    Cpx T_ref;
    double phase_ref, lambda_w;
    int i;

    build_path(fname, sizeof(fname), "fig3_phase_shift.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(5e3, 100e3, N, E2_arr);
    linspace(2.0, 50.0, N, eta2_arr);
    linspace(20e-6, 5e-6, N, h2_arr);

    T_ref = compute_T(f_center, h2_arr[0], E2_arr[0], eta2_arr[0]);
    phase_ref = cpx_arg(T_ref);
    lambda_w = c_w / f_center;

    fprintf(fp, "progress,E2_kPa,eta2_Pas,h2_um,T_abs,phase_deg,node_shift_um\n");
    for (i = 0; i < N; i++) {
        Cpx T = compute_T(f_center, h2_arr[i], E2_arr[i], eta2_arr[i]);
        double phase = cpx_arg(T) - phase_ref;
        double node_shift;
        while (phase > M_PI) phase -= 2*M_PI;
        while (phase < -M_PI) phase += 2*M_PI;
        node_shift = (phase / (2*M_PI)) * (lambda_w/2) * 1e6;
        fprintf(fp, "%.6f,%.3f,%.3f,%.3f,%.8f,%.4f,%.4f\n",
                (double)i/(N-1), E2_arr[i]/1e3, eta2_arr[i],
                h2_arr[i]*1e6, cpx_abs(T), phase*180.0/M_PI, node_shift);
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 4: Degradation evolution
 * ============================================================ */

static void output_degradation_evolution(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 1500, Nt = 6;
    double freqs[1500];
    double E2_t[6], eta2_t[6], h2_t[6];
    int i, t;

    build_path(fname, sizeof(fname), "fig4_degradation.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(17e6, 23e6, Nf, freqs);
    linspace(5e3, 50e3, Nt, E2_t);
    linspace(2.0, 20.0, Nt, eta2_t);
    linspace(20e-6, 10e-6, Nt, h2_t);

    fprintf(fp, "freq_MHz");
    for (t = 0; t < Nt; t++) fprintf(fp, ",T_t%d", t);
    fprintf(fp, "\n");

    for (i = 0; i < Nf; i++) {
        fprintf(fp, "%.6f", freqs[i]/1e6);
        for (t = 0; t < Nt; t++) {
            Cpx T = compute_T(freqs[i], h2_t[t], E2_t[t], eta2_t[t]);
            fprintf(fp, ",%.8f", cpx_abs(T));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 5: Attenuation and phase velocity
 * ============================================================ */

static void output_attenuation(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 500;
    double freqs[500];
    double eta_vals[] = {0.5, 2.0, 10.0, 50.0};
    int Ne = 4;
    int i, j;

    build_path(fname, sizeof(fname), "fig5_attenuation.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1e6, 50e6, Nf, freqs);

    fprintf(fp, "freq_MHz");
    for (j = 0; j < Ne; j++)
        fprintf(fp, ",alpha_dBmm_eta%.1f,c_phase_eta%.1f", eta_vals[j], eta_vals[j]);
    fprintf(fp, "\n");

    for (i = 0; i < Nf; i++) {
        double omega = 2.0 * M_PI * freqs[i];
        fprintf(fp, "%.6f", freqs[i]/1e6);
        for (j = 0; j < Ne; j++) {
            Cpx c_star = kv_wave_speed(omega, E_2_0, eta_vals[j], rho_2);
            Cpx k_star = cpx_div(cpx(omega, 0), c_star);
            double alpha = k_star.im;
            double alpha_dBmm = alpha * 20.0/log(10.0) * 1e-3;
            double c_phase = omega / k_star.re;
            fprintf(fp, ",%.6f,%.4f", alpha_dBmm, c_phase);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 6: Broadband sweep (1-100 MHz)
 * ============================================================ */

static void output_broadband_sweep(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 4000;
    double freqs[4000];
    double h4_vals[] = {330.0e-6, 165.0e-6, 110.0e-6, 82.0e-6, 66.0e-6};
    double h2_coupl = 10e-6;
    int Nh = 5;
    int i, j;

    build_path(fname, sizeof(fname), "fig6_broadband_sweep.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1.0e6, 100.0e6, Nf, freqs);

    fprintf(fp, "freq_MHz");
    for (j = 0; j < Nh; j++) {
        double f_res = c_4 / (2.0 * h4_vals[j]);
        fprintf(fp, ",T_h4_%.0fum_fres_%.1fMHz", h4_vals[j]*1e6, f_res/1e6);
    }
    fprintf(fp, "\n");

    for (i = 0; i < Nf; i++) {
        fprintf(fp, "%.6f", freqs[i]/1e6);
        for (j = 0; j < Nh; j++) {
            Cpx T = compute_T_full(freqs[i], h2_coupl, E_2_0, eta_2_0, h_3, h4_vals[j]);
            fprintf(fp, ",%.8f", cpx_abs(T));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 7: Optimal couplant thickness vs frequency
 * ============================================================ */

static void output_optimal_couplant(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 500;
    double freqs[500];
    int Nh2 = 1000;
    double h2_arr[1000];
    int i, j;

    build_path(fname, sizeof(fname), "fig7_optimal_couplant.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1.0e6, 80.0e6, Nf, freqs);
    linspace(0.1e-6, 200e-6, Nh2, h2_arr);

    fprintf(fp, "freq_MHz,h2_opt_um,h2_opt_over_lambda,T_max,"
                "T_at_1um,T_at_5um,T_at_20um,T_at_50um\n");

    for (i = 0; i < Nf; i++) {
        double T_max = 0.0;
        double h2_best = 0.0;
        double lambda2 = c_2 / freqs[i];
        double T_1, T_5, T_20, T_50;

        for (j = 0; j < Nh2; j++) {
            Cpx T = compute_T(freqs[i], h2_arr[j], E_2_0, eta_2_0);
            double Tabs = cpx_abs(T);
            if (Tabs > T_max) { T_max = Tabs; h2_best = h2_arr[j]; }
        }

        T_1  = cpx_abs(compute_T(freqs[i], 1e-6,  E_2_0, eta_2_0));
        T_5  = cpx_abs(compute_T(freqs[i], 5e-6,  E_2_0, eta_2_0));
        T_20 = cpx_abs(compute_T(freqs[i], 20e-6, E_2_0, eta_2_0));
        T_50 = cpx_abs(compute_T(freqs[i], 50e-6, E_2_0, eta_2_0));

        fprintf(fp, "%.6f,%.4f,%.6f,%.8f,%.8f,%.8f,%.8f,%.8f\n",
                freqs[i]/1e6, h2_best*1e6, h2_best/lambda2,
                T_max, T_1, T_5, T_20, T_50);
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 8: Glass resonance interaction map
 * ============================================================ */

static void output_glass_resonance_map(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 500, Nh3 = 200;
    double freqs[500], h3_arr[200];
    double h2_fixed = 10e-6;
    int i, j;

    build_path(fname, sizeof(fname), "fig8_glass_resonance.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(5.0e6, 60.0e6, Nf, freqs);
    linspace(0.1e-3, 3.0e-3, Nh3, h3_arr);

    fprintf(fp, "freq_MHz,h3_mm,T_abs,T_dB,glass_fres1_MHz\n");

    for (j = 0; j < Nh3; j++) {
        double fres1 = c_3 / (2.0 * h3_arr[j]);
        for (i = 0; i < Nf; i++) {
            Cpx T = compute_T_full(freqs[i], h2_fixed, E_2_0, eta_2_0, h3_arr[j], h_4);
            double Tabs = cpx_abs(T);
            double TdB = 10.0 * log10(dmax(Tabs, 1e-15));
            fprintf(fp, "%.6f,%.4f,%.8f,%.4f,%.4f\n",
                    freqs[i]/1e6, h3_arr[j]*1e3, Tabs, TdB, fres1/1e6);
        }
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 9: Multi-parameter design optimization
 * ============================================================ */

static void output_design_optimization(void)
{
    char fname[256];
    FILE *fp;
    double f_targets[] = {741e3, 8e5, 9e5, 1e6, 2e6, 5e6, 10e6, 20e6};
    /* Please change these f_targets as you wish */
    int Ntarg = 8;
    int Nh2 = 50, Nh3 = 30, Nh4 = 40;
    double h2_arr[50], h3_arr[30], h4_arr[40];
    int it, i2, i3, i4;

    build_path(fname, sizeof(fname), "fig9_design_optimization.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(0.5e-6, 100e-6, Nh2, h2_arr);
    linspace(0.1e-3, 2.0e-3, Nh3, h3_arr);
    linspace(30e-6, 500e-6, Nh4, h4_arr);

    fprintf(fp, "f_target_MHz,h2_opt_um,h3_opt_mm,h4_opt_um,"
                "T_max,f_res_transducer_MHz,f_res_glass_MHz,"
                "h2_over_lambda2,h3_over_lambda3\n");

    for (it = 0; it < Ntarg; it++) {
        double freq = f_targets[it];
        double T_best = 0.0;
        double h2_best = 0, h3_best = 0, h4_best = 0;

        if (freq >= 1e6)
            printf("  Optimizing f = %.0f MHz ... ", freq/1e6);
        else
            printf("  Optimizing f = %.0f kHz ... ", freq/1e3);
        fflush(stdout);

        for (i4 = 0; i4 < Nh4; i4++)
            for (i3 = 0; i3 < Nh3; i3++)
                for (i2 = 0; i2 < Nh2; i2++) {
                    Cpx T = compute_T_full(freq, h2_arr[i2], E_2_0, eta_2_0,
                                            h3_arr[i3], h4_arr[i4]);
                    double Tabs = cpx_abs(T);
                    if (Tabs > T_best) {
                        T_best = Tabs;
                        h2_best = h2_arr[i2];
                        h3_best = h3_arr[i3];
                        h4_best = h4_arr[i4];
                    }
                }

        fprintf(fp, "%.1f,%.4f,%.4f,%.4f,%.8f,%.4f,%.4f,%.6f,%.6f\n",
                freq/1e6, h2_best*1e6, h3_best*1e3, h4_best*1e6,
                T_best, c_4/(2.0*h4_best)/1e6, c_3/(2.0*h3_best)/1e6,
                h2_best / (c_2/freq), h3_best / (c_3/freq));

        printf("|T|_max = %.4f at h2=%.1f um, h3=%.2f mm, h4=%.0f um\n",
               T_best, h2_best*1e6, h3_best*1e3, h4_best*1e6);
    }
    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * MAIN
 * ============================================================ */

int main(void)
{
    double f_res_transducer, f_res_glass;
    int mkdir_result;

    /* Create output directory */
    mkdir_result = MKDIR(OUTPUT_DIR);
    if (mkdir_result != 0 && errno != EEXIST) {
        fprintf(stderr, "Error: cannot create %s\n", OUTPUT_DIR);
        return 1;
    }

    Z_1 = rho_1 * c_1_long;
    Z_3 = rho_3 * c_3;
    Z_4 = rho_4 * c_4;
    Z_w = rho_w * c_w;

    f_res_transducer = c_4 / (2.0 * h_4);
    f_res_glass = c_3 / (2.0 * h_3);

    printf("=====================================================\n");
    printf(" BENG207_1: Multilayer Acoustic Transfer-Matrix Model\n");
    printf(" Cross-platform build (Linux/Mac/Windows)\n");
    printf("=====================================================\n");
    printf(" Z_LiNbO3 (source) = %.2f MRayl\n", Z_1/1e6);
    printf(" Z_glass           = %.2f MRayl\n", Z_3/1e6);
    printf(" Z_LiNbO3 (cover)  = %.2f MRayl\n", Z_4/1e6);
    printf(" Z_water           = %.2f MRayl\n", Z_w/1e6);
    printf("-----------------------------------------------------\n");
    printf(" Transducer h4=%.0f um -> f_res = %.2f MHz\n",
           h_4*1e6, f_res_transducer/1e6);
    printf(" Glass      h3=%.0f mm -> f_res = %.2f MHz (n=1)\n",
           h_3*1e3, f_res_glass/1e6);
    printf(" lambda_w at 20 MHz   = %.1f um\n", c_w/f_center*1e6);
    printf(" lambda_coupl at 20 MHz = %.1f um\n", c_2/f_center*1e6);
    printf("=====================================================\n\n");

    printf("Generating original outputs (1-5)...\n\n");
    output_spectrum();
    output_sensitivity();
    output_phase_shift();
    output_degradation_evolution();
    output_attenuation();

    printf("\nGenerating broadband/optimization outputs (6-9)...\n\n");
    output_broadband_sweep();
    output_optimal_couplant();
    output_glass_resonance_map();
    output_design_optimization();

    printf("\nDone. 9 CSV files generated in %s/\n", OUTPUT_DIR);

#ifdef _WIN32
    printf("Press Enter to exit...");
    getchar();
#endif

    return 0;
}
