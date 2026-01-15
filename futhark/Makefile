# Define your ranges here
<<<<<<< HEAD
NUM_QUBITS := 5 200 400 800 1600 3200 6400 12000 24000 48000 96000 192000 
NUM_GATES  := 10 1000 10000 100000 1000000
=======
NUM_QUBITS := 200 400 800 1600 3200 6400 12000 24000 48000 96000 
NUM_GATES  := 10 100 1000 10000 100000 
>>>>>>> main

.PHONY: datasets
datasets:
	mkdir -p data
	@for q in $(NUM_QUBITS); do \
		midpoint=$$(($$q / 2)); \
		for g in $(NUM_GATES); do \
			echo "Generating data/$${q}_$${g}.in ..."; \
<<<<<<< HEAD
			python3 -c "import random, sys; \
			q = int(sys.argv[1]); \
			g = int(sys.argv[2]); \
			print(random.randint(0, 999999)); \
			print(q); \
			print(str([random.randint(0, 3) for _ in range(g)])); \
			q1s = [random.randint(0, q-1) for _ in range(g)]; \
			print(str(q1s)); \
			q2s = [(x + random.randint(1, q-1)) % q for x in q1s]; \
			print(str(q2s))" $$q $$g > data/$${q}_$${g}.in; \
=======
			futhark dataset \
				-g i32 \
				--i64-bounds=$${q}:$${q} -g i64 \
				--i64-bounds=0:3 -g [$${g}]i64 \
				--i64-bounds=0:$${midpoint} -g [$${g}]i64 \
				--i64-bounds=$${midpoint}:$${q} -g [$${g}]i64 \
				> data/$${q}_$${g}.in; \
>>>>>>> main
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
