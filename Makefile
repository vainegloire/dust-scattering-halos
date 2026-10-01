# Rebuild everything from the code.  Run from the repository root.
#
#   make            all figures and paper.pdf, from the stored results  (~2 min)
#   make test       unit test for the photon yields + ring-geometry checks
#   make figures    regenerate the eleven figures from results/
#   make paper      rebuild paper.tex and paper.pdf (needs pandoc + pdflatex)
#   make results    rerun every Monte Carlo from scratch  (~10 min)
#
# Every figure script fixes its random seed, so `make figures` reproduces the
# committed figures exactly from the stored results.  `make results`
# OVERWRITES results/ (compare with `git diff results/` afterwards).  It
# reproduces exactly: the multi-screen runs (the seeds of the 28 September 2026
# regeneration), the detection fractions (Fig. 9) and the calibration run.
# The original seeds of Table 1's single-screen columns, the single-screen
# self-calibration and the background sweep (Fig. 4) were not recorded, so
# those reproduce within Monte Carlo noise.  matched_filter.py reuses its
# seeded cache, results/matched_filter_cache.npz; delete that file to
# recompute it too (slow).

PY ?= python3
FIG_SCRIPTS = make_fig1_2.py crb_feasibility.py plot_fig4.py make_new_figs.py \
              fig7_realistic.py fig8_time_evolution.py fig9_plot.py matched_filter.py

.PHONY: all figures paper test results

all: figures paper

figures:
	cd code && for s in $(FIG_SCRIPTS); do echo "== $$s"; $(PY) $$s || exit 1; done

paper:
	bash build.sh

test:
	cd code && $(PY) test_yields.py && $(PY) sanity_check.py

results:
	cd results && rm -f mc_acc_ideal.json mc_acc_bkg.json mc_acc_multi.json \
	    selfcal_acc_single.json selfcal_acc_multi.json sightline_ensemble.json \
	    vsbkg_acc.json detect.json
	cd code && \
	  $(PY) mc_accumulate.py run ideal 120 2801 && \
	  $(PY) mc_accumulate.py run bkg 120 2802 && \
	  $(PY) mc_accumulate.py run multi 120 2809 && \
	  $(PY) mc_accumulate.py finalize && \
	  $(PY) selfcal_mc.py run single 120 2803 && \
	  $(PY) selfcal_mc.py run multi 120 2810 && \
	  $(PY) selfcal_mc.py finalize && \
	  $(PY) ensemble_sightlines.py 50 2811 && \
	  $(PY) ensemble_sightlines.py report && \
	  $(PY) vsbkg_smooth.py 2804 135 && \
	  $(PY) detect_mc.py 0.07 && \
	  $(PY) detect_mc.py 0.3 && \
	  $(PY) calib.py
