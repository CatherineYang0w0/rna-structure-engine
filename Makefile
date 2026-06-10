IMAGE ?= rna-structure-engine:stage0
LOCAL_IMAGE ?= uca1-rna-rnastructure:latest
CASE ?= xist
TEST ?= test-1
SCHEME ?= rnastr-deigan
REAGENT ?= 1M7
WASHIETL_SAMPLE_SIZE ?= 5
RUN ?= docker run --rm -v $(CURDIR):/work -w /work $(IMAGE)
RUN_LOCAL ?= docker run --rm -v $(CURDIR):/work -w /work $(LOCAL_IMAGE)

.PHONY: build run run-local stage1 stage1-local shell clean

build:
	docker build -t $(IMAGE) .

run:
	$(RUN) python -m engine.run_rnastructure \
		--case $(CASE) \
		--test $(TEST) \
		--scheme $(SCHEME) \
		--reagent $(REAGENT)

run-local:
	$(RUN_LOCAL) python -m engine.run_rnastructure \
		--case $(CASE) \
		--test $(TEST) \
		--scheme $(SCHEME) \
		--reagent $(REAGENT)

stage1:
	$(RUN) python -m engine.run_stage1 \
		--case $(CASE) \
		--test $(TEST) \
		--reagent $(REAGENT) \
		--washietl-sample-size $(WASHIETL_SAMPLE_SIZE)

stage1-local:
	$(RUN_LOCAL) python -m engine.run_stage1 \
		--case $(CASE) \
		--test $(TEST) \
		--reagent $(REAGENT) \
		--washietl-sample-size $(WASHIETL_SAMPLE_SIZE)

shell:
	$(RUN) bash

clean:
	rm -rf outputs/$(TEST)/$(CASE)
