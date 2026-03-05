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
	g++ -O2 -std=c++11 -I./include test/joye_libert_test.cpp src/*.cpp -o $(OUTDIR)/jltest -pthread -lntl -lgmp -lm

rinocchio: | out
	g++ -O2 -std=c++11 -I./include test/rinocchio_test.cpp src/*.cpp -o $(OUTDIR)/rinocchio -pthread -lntl -lgmp -lm

invertible: | out
	g++ -O2 -std=c++11 -I./include test/invertible_test.cpp src/*.cpp -o $(OUTDIR)/invertible -pthread -lntl -lgmp -lm

interpolate: | out
	g++ -O2 -std=c++11 -I./include test/interpolate_test.cpp src/*.cpp -o $(OUTDIR)/interpolate -pthread -lntl -lgmp -lm

out:
	mkdir -p $(OUTDIR)

# ============================================================
# Brakedown over GR targets
# ============================================================
BRAKEDOWN_SRCS = src/gr.cpp src/sparse_matrix_gr.cpp src/brakedown_code_gr.cpp \
                 src/merkle.cpp src/brakedown_pcs_gr.cpp

brakedown_test: | out
	g++ -O2 -std=c++11 -I./include test/brakedown_gr_test.cpp $(BRAKEDOWN_SRCS) \
		-o $(OUTDIR)/brakedown_gr_test -pthread -lntl -lgmp -lm

runbrakedown: brakedown_test
	./$(OUTDIR)/brakedown_gr_test

# ============================================================
# Fair Benchmark: 在同一进程中运行 Rinocchio + Brakedown
# 使用矩阵乘法电路 (与原始 Rinocchio 一致)
# ============================================================
fair_benchmark: | out
	g++ -O2 -std=c++11 -I./include test/fair_benchmark.cpp src/*.cpp \
		-o $(OUTDIR)/fair_benchmark -pthread -lntl -lgmp -lm

runfair: fair_benchmark
	./$(OUTDIR)/fair_benchmark

# ============================================================
# Rinocchio-only benchmark (独立运行)
# ============================================================
bench_commitment: | out
	g++ -O2 -std=c++11 -I./include test/bench_commitment.cpp src/*.cpp \
		-o $(OUTDIR)/bench_commitment -pthread -lntl -lgmp -lm

runbench: bench_commitment
	./$(OUTDIR)/bench_commitment