# Define your ranges here
NUM_QUBITS := 200 400 800 1600 3200 6400 12000 24000 48000 96000 
NUM_GATES  := 1000 10000 100000 1000000

.PHONY: datasets
datasets:
	mkdir -p data
	@for q in $(NUM_QUBITS); do \
		midpoint=$$(($$q / 2)); \
		for g in $(NUM_GATES); do \
			echo "Generating data/$${q}_$${g}.in ..."; \
			futhark dataset \
				-g i32 \  # random seed
				--i64-bounds=$${q}:$${q} -g i64 \ # num_qubits
				--i64-bounds=0:3 -g [$${g}]i64 \ # gate number
				--i64-bounds=0:$${midpoint} -g [$${g}]i64 \ # only allow control qubit to be lower half of possible qubits
				--i64-bounds=$${midpoint}:$${q} -g [$${g}]i64 \ # only allow target to be other half of qubits
				> data/$${q}_$${g}.in; \ # Output to data file
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
