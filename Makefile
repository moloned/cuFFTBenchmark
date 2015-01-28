ifeq ($(INPLACE), 1)
   USERFLAGS += -DINPLACE
endif
test: bench in.txt
	./bench < in.txt

bench: bench.cu
	nvcc -o bench bench.cu -DFROMTO=${FROMTO} ${USERFLAGS} -arch=sm_35 -lcufft
