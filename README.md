<pre> 
Title: BENG207_Multilayer_model
Author: BENG207 students (Instructor: Akihiro J. Matsuoka, M.D., D.M.Sc., Ph.D., FACS, Co-instructor: James Friend, Ph.D.)
Course: BENG207 Winter/Spring Semester (2026)
Status: Draft
Type: Analytical model
License: Public domain
Further discussion: TBD
</pre>

Instruction to BENG207 students
Please review all of the variables/parameters.

---------------

* [Abstract](#Abstract)
* [Background and Significance](#shopping-list)
* [Software Installation](#software-installation)
  * [Verifying your download](#verifying-your-download)
* [Enclosure Designs](#enclosure-designs)
* [SeedQR Printable Templates](#seedqr-printable-templates)
* [Build from Source](#build-from-source)
* [Developer Local Build Instructions](#developer-local-build-instructions)


---------------

## Abstract 
TBW

## Background and Significance
Our most recent iteration of the biohybrid cochlear implant, designed by A.J.M, resulted in fabrication of multiple prototype devices for feasibility testing. These experiments revealed a critical limitation in our prior approach to generating a static brain-derived neurotrophic factor (BDNF) concentration gradient. Although the initial gradient was sufficient to support neurite extension from human iPSC-derived otic neuronal progenitors and promote connectivity toward endogenous spiral ganglion neurons, the gradient was not temporally stable. By approximately day 7 in culture, the BDNF distribution became uniform within the microchannel, indicating dissipation of the intended concentration profile (Nella, et al., 2024).

The initial studies were conducted using a commercially available microfluidic platform (Zona), which was not designed to replicate the structural or biochemical microenvironment of the inner ear. Nevertheless, these experiments established proof-of-concept feasibility and provided a quantitative foundation for further device optimization.

Building upon these findings, we then initiated a second-generation design effort in collaboration with the BENG207 (2025) student cohort between January and May 2025. The revised microfluidic architecture was engineered to better approximate key aspects of the inner ear environment, including spatial confinement, gradient stability, and controlled mass transport dynamics. This redesign aims to support sustained neurotrophic signaling and improved structural integration in subsequent biohybrid implant prototypes. The initial results are encourging and Keristen Russ is now finalizing the experiments for the publication. 

The ratialnlae of the third gneneration of our microfluidic device was deisgned bvy AJM and James Friend, Ph.D. (Previously at UCSD now the chair of Wahshington University St. Louis) in Fall 2025 in that we now incoporate surface acoustic wave technology to generate long-lasting BDNF concnetration gradient. Based on Jia et al., 2025, The transfer-matrix (or transmission-line) method for multilayer acoustics is well-established. It seems that 
swapping in a complex wave speed from the Kelvin-Voigt model is straightforward. Frequency-dependent attenuation α(ω) ∝ ω²η/(2ρc³) at leading order, which means their 50 MHz system will be far more sensitive to couplant viscoelasticity than lower-frequency devices — that's I think is a testable condtionn. According to Alexi, we are using way lower frequency (somewhere around 1 MHz, please ask him). I can think of three failure modes I can list (increased attenuation, phase drift, impedance mismatch change) are genuinely distinct and could in principle be disentangled experimentally by tracking node positions versus overall force amplitude over time.


## 

## Output FIgures: Captions
Fig 1 — Transmission spectrum |T(f)| from 15–25 MHz for four couplant thicknesses (5, 15, 30, 50 µm), fresh vs degraded US gel

Fig 2a/2b — Sensitivity maps |T| vs (h₂, η₂) at fixed E₂, and |T| vs (h₂, E₂) at fixed η₂ — two heatmaps showing which couplant parameter matters more or less

Fig 3 — Phase shift during degradation — tracks how US gel aging (stiffening + thinning) shifts the transmission phase and displaces the standing wave nodes

Fig 4 — Degradation evolution — six time snapshots of the spectrum (17–23 MHz, this is way too high sorry) as gel degrades from fresh to old

Fig 5 — Couplant attenuation — frequency-dependent attenuation (dB/mm) and phase velocity in the US gel for four viscosity values (0.5, 2, 10, 50 Pa·s)

Fig 6 — Broadband sweep |T(f)| from 1–100 MHz for five transducer thicknesses (330, 165, 110, 82, 66 µm), each targeting a different center frequency (10, 20, 30, 40, 50 MHz) — shows how 
changing h₄ opens up different frequency bands

Fig 7 — Optimal couplant thickness vs frequency — for each frequency from 1–80 MHz, sweeps h₂ and finds the thickness that maximizes |T|, also reports |T| at four fixed thicknesses (1, 5, 20, 50 µm) for comparison

Fig 8 — Glass resonance interaction map — 2D heatmap of |T(f, h₃)| showing how glass thickness (0.1–3 mm) creates periodic transmission windows at its resonance harmonics f_n = n·c₃/(2h₃)

Fig 9 — Multi-parameter design optimization — brute-force g;(AKA primitive) rid search over (h₂, h₃, h₄) simultaneously for eight target frequencies, reports optimal geometry and maximum |T| for each
