# WfAverage analysis workflow
# Usage:
#   make [all]                        cmake build (default)
#   make RUN=<id> build_list           generate EOS file list
#   make RUN=<id> analysis             run waveform averaging
#   make RUN=<id> MONITOR_CHAN=<n> analysis   (override monitor PMT channel)
#   make RUN=<id> plot                 draw averaged waveforms (±1σ band)
#   make RUN=<id> plot_no_band         draw averaged waveforms (no band)
#   make clean                         remove build/ and analysis/

OUT       ?= analysis
ANA_OPT   ?= --time-align --trigger-type="Calibration"
MONITOR_CHAN ?= 43303
EOS_BASE  ?= /eos/juno/juno-rtraw
JUNO_SW   ?= J25.7.1
TRIGGER   ?= global_trigger
EVTMAX    ?= -1

BUILD_DIR    = build
LIB_WFA      = $(BUILD_DIR)/lib/libWfAverage.so
DRAW_EXE     = $(BUILD_DIR)/bin/draw_wf_avg.exe
LIST_FILE    = $(OUT)/$(RUN).list
ROOT_FILE    = $(OUT)/wf_avg_$(RUN).root
PLOT_FILE    = $(OUT)/wf_avg_$(RUN)_plots.pdf
PLOT_NOBAND  = $(OUT)/wf_avg_$(RUN)_noband.pdf

.PHONY: all build build_list analysis plot plot_no_band clean check_run

all: build

build:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_EXPORT_COMPILE_COMMANDS=true
	cmake --build $(BUILD_DIR)

check_run:
	@if [ -z "$(RUN)" ]; then \
		echo "ERROR: RUN is not set. Usage: make RUN=<run_id> <target>"; \
		exit 1; \
	fi

# ---- EOS list generation ----
$(LIST_FILE): gen_run_list.sh | check_run
	@mkdir -p $(OUT)
	@if [ -n "$(EOS_PATTERN)" ]; then \
		pat="$(EOS_PATTERN)"; \
	else \
		sub=$$(printf '%08d' $$(( $(RUN) / 100 * 100 ))); \
		pat="$(EOS_BASE)/$(JUNO_SW)/$(TRIGGER)/$$sub/$$sub/$(RUN)/*.rtraw"; \
	fi; \
	./gen_run_list.sh "$$pat"
	@mv $(RUN).list $@
	@echo "List  -> $@"

build_list: check_run $(LIST_FILE)

# ---- Waveform averaging ----
$(ROOT_FILE): $(LIST_FILE) run.py WfAverage.cxx $(LIB_WFA)
	python run.py --input-list $(LIST_FILE) --evtmax $(EVTMAX) \
		--user-output $(ROOT_FILE) $(ANA_OPT) --monitor-channel $(MONITOR_CHAN)

analysis: check_run $(ROOT_FILE)

# ---- Plotting (with ±1σ band) ----
$(PLOT_FILE): $(ROOT_FILE) $(DRAW_EXE)
	$(DRAW_EXE) --input $(ROOT_FILE) --output $(PLOT_FILE)

plot: check_run $(PLOT_FILE)

# ---- Plotting (no band) ----
$(PLOT_NOBAND): $(ROOT_FILE) $(DRAW_EXE)
	$(DRAW_EXE) --input $(ROOT_FILE) --output $(PLOT_NOBAND) --no-band

plot_no_band: check_run $(PLOT_NOBAND)

clean:
	rm -rf $(BUILD_DIR) $(OUT)
