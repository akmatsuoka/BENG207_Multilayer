/* main.c
 * This is the least difficult model
 * Multilayer Acoustic Transfer-Matrix Model
 * with Kelvin-Voigt viscoelastic couplant (ultrasound gels)
 *
 * Stack (bottom to top) See also my Figure 1:
 *   1. LiNbO3 transducer (semi-infinite source)
 *   2. Ultrasound gel couplant (Kelvin-Voigt: E* = E + iωη)
 *   3. Borosilicate glass superstrate AKA a slide glass (e.g., 1 mm)
 *   4. LiNbO3 coverslip AKA a cover glass (e.g., 150 um)
 *   5. Water fluid channel (terminating half-space)
 *
 * Output: CSV files for plotting use whatever Anna likes to use.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Alias for imaginary unit to avoid conflict with variable names */
#define IMAG_UNIT (_Complex_I)
/* if you have a question, please let me know. I almost always do this in C to avoid pain.
 *
/* Output directory — change YOUR_USERNAME to your actual username */
#define OUTPUT_DIR "/home/shh/Projects/BENG207_1/PLOT_CSV/"

/* ============================================================
 * Material parameters
 * ============================================================ */

/* Layer 1: LiNbO3 transducer (source) */
static const double rho_1     = 4647.0;    /* kg/m^3 */
static const double c_1_long  = 6570.0;    /* m/s longitudinal */

/* Layer 2: US gel couplant (Kelvin-Voigt) */
static const double rho_2   = 1020.0;      /* kg/m^3 */
static const double E_2_0   = 5.0e3;       /* Pa (fresh gel) */
static const double eta_2_0 = 2.0;         /* Pa-s (fresh gel) */

/* Layer 3: Borosilicate glass */
static const double rho_3 = 2230.0;        /* kg/m^3 */
static const double c_3   = 5640.0;        /* m/s */
static const double h_3   = 1.0e-3;        /* 1 mm */

/* Layer 4: LiNbO3 coverslip */
static const double rho_4 = 4647.0;        /* kg/m^3 */
static const double c_4   = 6570.0;        /* m/s */
static const double h_4   = 150.0e-6;      /* 150 um */

/* Layer 5: Water (load) */
static const double rho_w = 1000.0;        /* kg/m^3 */
static const double c_w   = 1483.0;        /* m/s */

/* Actuation */
static const double f_center = 20.0e6;     /* 20 MHz */

/* Derived (set in main) */
static double Z_1, Z_3, Z_4, Z_w;

/* ============================================================
 * 2x2 complex matrix type
 * ============================================================ */

typedef struct {
    double complex m[2][2];
} Mat2x2;

/* Matrix multiply: C = A * B */
static Mat2x2 mat_mul(Mat2x2 A, Mat2x2 B)
{
    Mat2x2 C;
    int i, j;
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            C.m[i][j] = A.m[i][0]*B.m[0][j] + A.m[i][1]*B.m[1][j];
    return C;
}

/* ============================================================
 * Kelvin-Voigt complex wave speed
 * c*(omega) = sqrt( (E + i*omega*eta) / rho )
 * ============================================================ */

static double complex kv_wave_speed(double omega, double E, double eta,
                                     double rho)
{
    double complex E_star = E + IMAG_UNIT * omega * eta;
    return csqrt(E_star / rho);
}

/* ============================================================
 * Transfer matrix for a single layer
 *
 *       [ cos(kh)        i sin(kh)/Z ]
 * M  =  [                             ]
 *       [ i Z sin(kh)    cos(kh)      ]
 * ============================================================ */

static Mat2x2 layer_matrix(double omega, double rho, double complex c_cmplx,
                            double h)
{
    double complex k  = omega / c_cmplx;
    double complex Z  = rho * c_cmplx;
    double complex kh = k * h;
    double complex cos_kh = ccos(kh);
    double complex sin_kh = csin(kh);
    Mat2x2 M;

    M.m[0][0] = cos_kh;
    M.m[0][1] = IMAG_UNIT * sin_kh / Z;
    M.m[1][0] = IMAG_UNIT * Z * sin_kh;
    M.m[1][1] = cos_kh;
    return M;
}

/* ============================================================
 * Full stack transmission coefficient
 *
 * T = 2 / (M11 + M12*Z_load + M21/Z_src + M22*Z_load/Z_src)
 * ============================================================ */

static double complex compute_T(double freq, double h2, double E2,
                                 double eta2)
{
    double omega = 2.0 * M_PI * freq;
    double complex c2, denom;
    Mat2x2 M2, M3, M4, Mt;

    /* Layer 2: couplant (KV) */
    c2 = kv_wave_speed(omega, E2, eta2, rho_2);
    M2 = layer_matrix(omega, rho_2, c2, h2);

    /* Layer 3: glass */
    M3 = layer_matrix(omega, rho_3, c_3, h_3);

    /* Layer 4: coverslip */
    M4 = layer_matrix(omega, rho_4, c_4, h_4);

    /* Total: M = M2 * M3 * M4 */
    Mt = mat_mul(M2, M3);
    Mt = mat_mul(Mt, M4);

    denom = Mt.m[0][0]
          + Mt.m[0][1] * Z_w
          + Mt.m[1][0] / Z_1
          + Mt.m[1][1] * Z_w / Z_1;

    return 2.0 / denom;
}

/* ============================================================
 * Helper: linspace
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

/* ============================================================
 * OUTPUT 1: Transmission spectrum |T(f)| for various h2
 * ============================================================ */

static void output_spectrum(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 2000;
    double freqs[2000];
    double h2_vals[] = {5e-6, 15e-6, 30e-6, 50e-6};
    int Nh = 4;
    double E_deg  = 50.0e3;
    double eta_deg = 20.0;
    int i, j;

    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "fig1_spectrum.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(15.0e6, 25.0e6, Nf, freqs);

    /* Header */
    fprintf(fp, "freq_MHz");
    for (j = 0; j < Nh; j++)
        fprintf(fp, ",T_fresh_h%.0fum,T_deg_h%.0fum",
                h2_vals[j]*1e6, h2_vals[j]*1e6);
    fprintf(fp, "\n");

    for (i = 0; i < Nf; i++) {
        fprintf(fp, "%.6f", freqs[i] / 1e6);
        for (j = 0; j < Nh; j++) {
            double complex Tf = compute_T(freqs[i], h2_vals[j],
                                          E_2_0, eta_2_0);
            double complex Td = compute_T(freqs[i], h2_vals[j],
                                          E_deg, eta_deg);
            fprintf(fp, ",%.8f,%.8f", cabs(Tf), cabs(Td));
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("  Written: %s (%d freq x %d thicknesses)\n", fname, Nf, Nh);
}

/* ============================================================
 * OUTPUT 2: Sensitivity maps |T| vs (h2, eta2) and (h2, E2)
 * ============================================================ */

static void output_sensitivity(void)
{
    char fname1[256], fname2[256];
    int Nh = 100;
    int Np = 100;
    double h2_arr[100], eta_arr[100], E_arr[100];
    FILE *fp1, *fp2;
    int i, j;

    snprintf(fname1, sizeof(fname1), "%s%s", OUTPUT_DIR, "fig2a_sensitivity_eta.csv");
    snprintf(fname2, sizeof(fname2), "%s%s", OUTPUT_DIR, "fig2b_sensitivity_E.csv");

    linspace(1e-6, 50e-6, Nh, h2_arr);
    logspace(-1, 2, Np, eta_arr);    /* 0.1 to 100 Pa-s */
    logspace(2, 6, Np, E_arr);       /* 100 Pa to 1 MPa */

    /* Map 1: |T| vs (h2, eta2) at fixed E2 */
    fp1 = fopen(fname1, "w");
    if (!fp1) { perror(fname1); return; }
    fprintf(fp1, "h2_um,eta2_Pas,T_abs\n");
    for (i = 0; i < Np; i++) {
        for (j = 0; j < Nh; j++) {
            double complex T = compute_T(f_center, h2_arr[j],
                                         E_2_0, eta_arr[i]);
            fprintf(fp1, "%.4f,%.6f,%.8f\n",
                    h2_arr[j]*1e6, eta_arr[i], cabs(T));
        }
    }
    fclose(fp1);
    printf("  Written: %s (%dx%d)\n", fname1, Np, Nh);

    /* Map 2: |T| vs (h2, E2) at fixed eta2 */
    fp2 = fopen(fname2, "w");
    if (!fp2) { perror(fname2); return; }
    fprintf(fp2, "h2_um,E2_Pa,T_abs\n");
    for (i = 0; i < Np; i++) {
        for (j = 0; j < Nh; j++) {
            double complex T = compute_T(f_center, h2_arr[j],
                                         E_arr[i], eta_2_0);
            fprintf(fp2, "%.4f,%.2f,%.8f\n",
                    h2_arr[j]*1e6, E_arr[i], cabs(T));
        }
    }
    fclose(fp2);
    printf("  Written: %s (%dx%d)\n", fname2, Np, Nh);
}

/* ============================================================
 * OUTPUT 3: Phase shift and node displacement during degradation
 * ============================================================ */

static void output_phase_shift(void)
{
    char fname[256];
    FILE *fp;
    int N = 200;
    double E2_arr[200], eta2_arr[200], h2_arr[200];
    double complex T_ref;
    double phase_ref, lambda_w;
    int i;

    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "fig3_phase_shift.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(5e3,   100e3, N, E2_arr);
    linspace(2.0,   50.0,  N, eta2_arr);
    linspace(20e-6, 5e-6,  N, h2_arr);

    T_ref = compute_T(f_center, h2_arr[0], E2_arr[0], eta2_arr[0]);
    phase_ref = carg(T_ref);
    lambda_w  = c_w / f_center;

    fprintf(fp, "progress,E2_kPa,eta2_Pas,h2_um,"
                "T_abs,phase_deg,node_shift_um\n");

    for (i = 0; i < N; i++) {
        double complex T = compute_T(f_center, h2_arr[i],
                                      E2_arr[i], eta2_arr[i]);
        double phase = carg(T) - phase_ref;
        double node_shift;

        /* Unwrap */
        while (phase >  M_PI) phase -= 2*M_PI;
        while (phase < -M_PI) phase += 2*M_PI;

        node_shift = (phase / (2*M_PI)) * (lambda_w/2) * 1e6;

        fprintf(fp, "%.6f,%.3f,%.3f,%.3f,%.8f,%.4f,%.4f\n",
                (double)i/(N-1),
                E2_arr[i]/1e3, eta2_arr[i], h2_arr[i]*1e6,
                cabs(T), phase*180.0/M_PI, node_shift);
    }

    fclose(fp);
    printf("  Written: %s (%d steps)\n", fname, N);
}

/* ============================================================
 * OUTPUT 4: Degradation time evolution of spectrum
 * ============================================================ */

static void output_degradation_evolution(void)
{
    char fname[256];
    FILE *fp;
    int Nf = 1500;
    int Nt = 6;
    double freqs[1500];
    double E2_t[6], eta2_t[6], h2_t[6];
    int i, t;

    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "fig4_degradation.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(17e6, 23e6, Nf, freqs);
    linspace(5e3,  50e3,  Nt, E2_t);
    linspace(2.0,  20.0,  Nt, eta2_t);
    linspace(20e-6, 10e-6, Nt, h2_t);

    /* Header */
    fprintf(fp, "freq_MHz");
    for (t = 0; t < Nt; t++)
        fprintf(fp, ",T_t%d", t);
    fprintf(fp, "\n");

    for (i = 0; i < Nf; i++) {
        fprintf(fp, "%.6f", freqs[i]/1e6);
        for (t = 0; t < Nt; t++) {
            double complex T = compute_T(freqs[i], h2_t[t],
                                          E2_t[t], eta2_t[t]);
            fprintf(fp, ",%.8f", cabs(T));
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("  Written: %s (%d freq x %d snapshots)\n", fname, Nf, Nt);
}

/* ============================================================
 * OUTPUT 5: Attenuation and phase velocity in couplant
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

    snprintf(fname, sizeof(fname), "%s%s", OUTPUT_DIR, "fig5_attenuation.csv");
    fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    linspace(1e6, 50e6, Nf, freqs);

    fprintf(fp, "freq_MHz");
    for (j = 0; j < Ne; j++)
        fprintf(fp, ",alpha_dBmm_eta%.1f,c_phase_eta%.1f",
                eta_vals[j], eta_vals[j]);
    fprintf(fp, "\n");

    for (i = 0; i < Nf; i++) {
        double omega = 2.0 * M_PI * freqs[i];
        fprintf(fp, "%.6f", freqs[i]/1e6);

        for (j = 0; j < Ne; j++) {
            double complex c_star = kv_wave_speed(omega, E_2_0,
                                                   eta_vals[j], rho_2);
            double complex k_star = omega / c_star;
            double alpha = cimag(k_star);
            double alpha_dBmm = alpha * 20.0/log(10.0) * 1e-3;
            double c_phase = omega / creal(k_star);

            fprintf(fp, ",%.6f,%.4f", alpha_dBmm, c_phase);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("  Written: %s (%d freq x %d viscosities)\n", fname, Nf, Ne);
}

/* ============================================================
 * MAIN
 * ============================================================ */

int main(void)
{
    /* Create output directory if it doesn't exist */
    if (mkdir(OUTPUT_DIR, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "Error: cannot create %s\n", OUTPUT_DIR);
        return 1;
    }

    /* Compute derived impedances */
    Z_1 = rho_1 * c_1_long;
    Z_3 = rho_3 * c_3;
    Z_4 = rho_4 * c_4;
    Z_w = rho_w * c_w;

    printf("=====================================================\n");
    printf(" Multilayer Acoustic Transfer-Matrix Model (KV)\n");
    printf("=====================================================\n");
    printf(" Z_LiNbO3 (source) = %.2f MRayl\n", Z_1/1e6);
    printf(" Z_glass           = %.2f MRayl\n", Z_3/1e6);
    printf(" Z_LiNbO3 (cover)  = %.2f MRayl\n", Z_4/1e6);
    printf(" Z_water           = %.2f MRayl\n", Z_w/1e6);
    printf(" lambda_w at 20MHz = %.1f um\n", c_w/f_center*1e6);
    printf(" Glass res spacing = %.2f MHz\n", c_3/(2*h_3)/1e6);
    printf("=====================================================\n\n");

    printf("Generating output files...\n\n");

    output_spectrum();
    output_sensitivity();
    output_phase_shift();
    output_degradation_evolution();
    output_attenuation();

    printf("\nDone. Plot CSVs with gnuplot, Python, or your tool of choice.\n");

    return 0;
}
