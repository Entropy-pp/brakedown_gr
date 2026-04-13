OUTDIR = out

build: joye_libert rinocchio

test: joye_libert
	./$(OUTDIR)/jltest

run: rinocchio
	./$(OUTDIR)/rinocchio

runinv: invertible
	./$(OUTDIR)/invertible

runint: interpolate
	./$(OUTDIR)/interpolate

joye_libert: | out
	g++ -O2 -std=c++11 -I./include test/joye_libert_test.cpp src/*.cpp -o $(OUTDIR)/jltest -pthread -lntl -lgmp -lm -lssl -lcrypto

rinocchio: | out
	g++ -O2 -std=c++11 -I./include test/rinocchio_test.cpp src/*.cpp -o $(OUTDIR)/rinocchio -pthread -lntl -lgmp -lm -lssl -lcrypto

invertible: | out
	g++ -O2 -std=c++11 -I./include test/invertible_test.cpp src/*.cpp -o $(OUTDIR)/invertible -pthread -lntl -lgmp -lm -lssl -lcrypto

interpolate: | out
	g++ -O2 -std=c++11 -I./include test/interpolate_test.cpp src/*.cpp -o $(OUTDIR)/interpolate -pthread -lntl -lgmp -lm -lssl -lcrypto

out:
	mkdir -p $(OUTDIR)

# ============================================================
# Brakedown over GR targets
# ============================================================
BRAKEDOWN_SRCS = src/gr.cpp src/sparse_matrix_gr.cpp src/brakedown_code_gr.cpp \
                 src/merkle.cpp src/brakedown_pcs_gr.cpp

brakedown_test: | out
	g++ -O2 -std=c++11 -I./include test/brakedown_gr_test.cpp $(BRAKEDOWN_SRCS) \
		-o $(OUTDIR)/brakedown_gr_test -pthread -lntl -lgmp -lm -lssl -lcrypto

runbrakedown: brakedown_test
	./$(OUTDIR)/brakedown_gr_test

# ============================================================
# Fair Benchmark
# ============================================================
fair_benchmark: | out
	g++ -O2 -std=c++11 -I./include test/fair_benchmark.cpp src/*.cpp \
		-o $(OUTDIR)/fair_benchmark -pthread -lntl -lgmp -lm -lssl -lcrypto

runfair: fair_benchmark
	./$(OUTDIR)/fair_benchmark

# ============================================================
# Rinocchio-only benchmark
# ============================================================
bench_commitment: | out
	g++ -O2 -std=c++11 -I./include test/bench_commitment.cpp src/*.cpp \
		-o $(OUTDIR)/bench_commitment -pthread -lntl -lgmp -lm -lssl -lcrypto

runbench: bench_commitment
	./$(OUTDIR)/bench_commitment

# ============================================================
# Brakedown PCS Benchmark
# ============================================================
bench_pcs: | out
	g++ -O2 -std=c++11 -I./include test/bench_pcs.cpp $(BRAKEDOWN_SRCS) \
		-o $(OUTDIR)/bench_pcs -pthread -lntl -lgmp -lm -lssl -lcrypto

runbenchpcs: bench_pcs
	./$(OUTDIR)/bench_pcs

# ============================================================
# Small Ring PCS Test
# ============================================================
small_ring_test: | out
	g++ -O2 -std=c++11 -I./include test/small_ring_test.cpp $(BRAKEDOWN_SRCS) \
		-o $(OUTDIR)/small_ring_test -pthread -lntl -lgmp -lm -lssl -lcrypto

runsmallring: small_ring_test
	./$(OUTDIR)/small_ring_test

# ============================================================
# Full Benchmark: Large Ring + Small Ring
# ============================================================
bench_brakedown: | out
	g++ -O2 -std=c++11 -I./include test/bench_brakedown.cpp $(BRAKEDOWN_SRCS) \
		-o $(OUTDIR)/bench_brakedown -pthread -lntl -lgmp -lm -lssl -lcrypto

runbenchbrakedown: bench_brakedown
	./$(OUTDIR)/bench_brakedown

# ============================================================
# Distributed Brakedown Benchmark
# ============================================================
DIST_SRCS = src/gr.cpp src/sparse_matrix_gr.cpp src/brakedown_code_gr.cpp \
            src/merkle.cpp src/brakedown_pcs_gr.cpp src/brakedown_distributed.cpp

bench_distributed: | out
	g++ -O3 -std=c++17 -pthread -I./include test/bench_distributed.cpp $(DIST_SRCS) \
		-o $(OUTDIR)/bench_distributed -lntl -lgmp -lm -lssl -lcrypto

rundist: bench_distributed
	./$(OUTDIR)/bench_distributed