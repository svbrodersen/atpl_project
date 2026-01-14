# Define your ranges here
NUM_QUBITS := 5 200 400 800 1600 3200 6400 12000 24000 48000 96000 192000 
NUM_GATES  := 10 1000 10000 100000 1000000

.PHONY: datasets
datasets:
	mkdir -p data
	@for q in $(NUM_QUBITS); do \
		for g in $(NUM_GATES); do \
			echo "Generating data/$${q}_$${g}.in ..."; \
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
		done \
	done
