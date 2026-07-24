//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef TENSOR_H_
#define TENSOR_H_

#include <iostream>
#include <memory>
#include <string>
#include <torch/script.h>

//----------------------------------------------------------
// LibTorch (torch::jit) inference session.
//
// DNAscent formerly used the TensorFlow C API to run a SavedModel.  The
// detection model is now a PyTorch Mamba model exported as a TorchScript
// module (model.pt) via utils/export_torchscript.py.  ModelSession wraps a
// torch::jit::script::Module together with the device the module lives on.
//
// The scripted module was saved in eval() mode, so its forward() uses a
// memory-efficient sequential scan that never materialises a (B, L, D, N)
// hidden-state tensor.  This makes inference safe for reads that are hundreds
// of kb long.
//----------------------------------------------------------
class ModelSession{
	public:
		torch::jit::script::Module module;
		torch::Device device = torch::Device(torch::kCPU);
};

// Load the TorchScript module on CPU.  'threads' sets the intra-op thread pool.
std::shared_ptr<ModelSession> model_load_cpu(const char *model_path, unsigned int threads);

// Load the TorchScript module on the given CUDA device.
std::shared_ptr<ModelSession> model_load_gpu(const char *model_path, unsigned char device, unsigned int threads);

#endif
