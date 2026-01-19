.PHONY: futhark backend
futhark:
	$(MAKE) -C futhark benchmark

backend:
	$(MAKE) -C futhark backend

hqp: backend
	cd interface/HQPlayground/; cabal build

# required backend targets for cabal interface builds
interface: hqp
	cd interface/HQPlayground/; cabal run interface

hybrid: hqp
	cd interface/HQPlayground/; cabal run HybridGPUCPU

unsupported: hqp
	cd interface/HQPlayground/; cabal run UnsupportedGPU

GPUHQP: hqp
	cd interface/HQPlayground/; cabal run GPUHQPTests

python: 
	pip install numpy pyperf qiskit qiskit-aer-gpu
	cd ..; python3 run_benchmarks.py

clean:
	rm -rf interface/HQPlayground/futhark-gen/hqp_gpu*
	cd interface/HQPlayground/; cabal clean
	$(MAKE) -C futhark clean
