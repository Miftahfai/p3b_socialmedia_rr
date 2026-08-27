# Social-media use and the occipital P3b in adolescents — Registered Report (Stage 1)

Analysis code for the pilot study and the power analysis supporting the Stage 1
Registered Report. **No participant-level data are included in this repository.**
The pilot data (adolescent participants) are held in secure institutional storage
at the University of Leeds and will be deposited in a public repository, with a
permanent archived version and DOI, upon acceptance of the Stage 2 manuscript.

## Contents

    power/
      rr_power_spearman.R                        Sample-size calculation
                                                 (safeguard anchor; May & Looney, 2020).
                                                 Runs standalone — no data required.
    preprocessing/
      Preprocessing_allch_p3b.ipynb              EEG preprocessing (MNE + autoreject)
      extract_behaviour_p3b.R                    behavioural accuracy / RT
      timing_rt_check.py                         stimulus-timing checks
      autoreject_interpolation_summary.py
      p3b_bad_channel_check.ipynb
    analysis/
      p3b_composite_O1O2avg.R                    Primary pilot analysis + pilot figures
      p3b_composite_erp_corrtopo_O1O2avg.ipynb   Grand-average ERP + P3b topography

## Requirements

**R (>= 4.4)** — tidyverse, MASS, performance, dagitty

    install.packages(c("tidyverse", "MASS", "performance", "dagitty"))

**Python (>= 3.10)** — mne, autoreject, numpy, pandas, scipy, scikit-learn, matplotlib

    pip install mne autoreject numpy pandas scipy scikit-learn matplotlib

## Reproduce

    Rscript power/rr_power_spearman.R      # power analysis — no data required

The preprocessing and pilot-analysis scripts require the participant-level data,
which are not shared at this stage (see Data availability in the manuscript).

## Data availability

Participant-level pilot data are held in secure institutional storage at the
University of Leeds and are available to editors and reviewers on request during
review. They will be deposited in a public repository and shared publicly upon
acceptance of the Stage 2 manuscript.
