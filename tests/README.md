# TensorFlow configuration regression test

This test needs a C++14 compiler and the TensorFlow C library already used by
DNAscent (`tensorflow/include` and `tensorflow/lib`). It does not load a model,
create compute threads, or require a GPU or the sequencing-data dependencies.

From the repository root:

```sh
mkdir -p build
g++ -std=c++14 -Wall -Wextra -Werror -I tensorflow/include \
    tests/tensor_config.cpp -L tensorflow/lib \
    -Wl,-rpath,"$PWD/tensorflow/lib" -ltensorflow -o build/test_tensor_config
build/test_tensor_config
```

The test checks golden protobuf messages, byte-for-byte compatibility for small
thread counts, and acceptance by `TF_SetConfig` for CPU and GPU configurations
across the one-, two-, and three-byte varint boundaries, up to `INT_MAX`.
GPU selectors `0`, `1`, and `9` are tested. A failed check exits nonzero; successful
cases print the thread count, ASCII device selector (zero for CPU), and serialized
ConfigProto in hexadecimal for independent decoding.
