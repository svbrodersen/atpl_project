# Define your ranges here
NUM_QUBITS := 5 200 400 800 1600 3200 6400 12000 24000 48000 96000 
NUM_GATES  := 10 100 1000 10000 100000 

.PHONY: datasets
datasets:
	mkdir -p data
	@for q in $(NUM_QUBITS); do \
		midpoint=$$(($$q / 2)); \
		max_val=$$(($$q - 1)); \
		for g in $(NUM_GATES); do \
			echo "Generating data/$${q}_$${g}.in ..."; \
			futhark dataset \
				-g i32 \
				--i64-bounds=$${q}:$${q} -g i64 \
				--i64-bounds=0:3 -g [$${g}]i64 \
				--i64-bounds=0:$${midpoint} -g [$${g}]i64 \
				--i64-bounds=$${midpoint}:$${max_val} -g [$${g}]i64 \
				> data/$${q}_$${g}.in; \
		done \
	done

benchmark:
	mkdir -p results
	futhark bench --backend=cuda --json=results/results.json main.fut > results/result.out

futhark_comment:
	@for n in $(NUM_QUBITS); do \
		for m in $(NUM_GATES); do \
				echo "-- notest compiled input @ data/$${n}_$${m}.in"; \
			done \
	done
