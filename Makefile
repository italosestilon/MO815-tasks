FLAGS = -fPIC -std=gnu11 -Wall -Wno-unused-result -pedantic

BIN = ./bin

LIBSVM_DIR = ./externals/libsvm
LIBNIFTI_DIR= ./externals/libnifti
LIBJPEG_DIR= ./externals/libjpeg
LIBPNG_DIR= ./externals/libpng
TSNE_DIR= ./externals/tsne

LIBIFT_INC   = -I ./include
LIBIFT_LD    = -L ./lib -lift
LIBSVM_INC   = -I $(LIBSVM_DIR)/include 
LIBSVM_LD    = -L $(LIBSVM_DIR)/lib -lsvm -lstdc++ 
LIBCBLAS_INC = -I /usr/local/opt/openblas/include
LIBCBLAS_LD  = -L /usr/local/opt/openblas/lib -L /lib64/atlas-sse3 -L /usr/lib/atlas-base -llapack -lblas -lcblas

LIBNIFTI_INC = -I $(LIBNIFTI_DIR)/include
LIBJPEG_INC = -I $(LIBJPEG_DIR)/include
TSNE_INC = -I $(TSNE_DIR)/include

EXTERNALS_LD = -fopenmp -lm -lpng -lz

INCLUDES = $(LIBIFT_INC) $(LIBSVM_INC) $(LIBCBLAS_INC) $(LIBNIFTI_INC) $(LIBJPEG_INC) $(TSNE_INC)
LIBS     = $(LIBIFT_LD) $(LIBSVM_LD) $(LIBCBLAS_LD) $(EXTERNALS_LD)

#IFT_DEBUG=1
ifeq ($(IFT_DEBUG), 1)
    export FLAGS += -pg -g -fsanitize=address -fsanitize=leak -DIFT_DEBUG=1
else
    export FLAGS += -O3
endif

#CUDA path in case IFT_GPU is enabled
export CUDA_DIR1=/usr/local/cuda
export CUDA_DIR2=/opt/cuda

ifeq ($(IFT_GPU), 1)
    export FLAGS += -DIFT_GPU=1
    INCLUDES += -I $(CUDA_DIR1)/include
    INCLUDES += -I $(CUDA_DIR2)/include
    LIBIFT_LD   += -L $(CUDA_DIR1)/lib64 -L $(CUDA_DIR2)/lib64 -lcublas -lcudart
endif



$@.c: $@.c
	$(CC) $(FLAGS) $@.c -o $(BIN)/$@ $(INCLUDES) $(LIBS)


clean:
	rm -rf ./bin/*










