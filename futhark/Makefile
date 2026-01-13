# Define your ranges here
NUM_QUBITS := 200 400 800 1600 3200 6400 12000 24000 48000 96000 192000 384000
NUM_GATES  := 1000 10000 100000 1000000

.PHONY: datasets
datasets:
	mkdir -p data
	@for q in $(NUM_QUBITS); do \
		for g in $(NUM_GATES); do \
			echo "Generating data/$${q}_$${g}.in ..."; \
			futhark dataset \
				-g i32 --i64-bounds=$${q}:$${q} \
				-g i64 --i64-bounds=0:3 \
				-g [$${g}]i64 --i64-bounds=0:$${q} \
				-g [$${g}]i64 --i64-bounds=0:$${q} \
				-g [$${g}]i64 \
				> data/$${q}_$${g}.in; \
		done \
	done
