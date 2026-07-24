//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include "tensor.h"
#include <torch/torch.h>
#include <torch/cuda.h>
#include <torch/csrc/jit/runtime/graph_executor.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>


static std::shared_ptr<ModelSession> load_module(const char *model_path, torch::Device device, unsigned int threads){

	// DNAscent parallelises at the READ level: `detect` launches `threads`
	// OpenMP workers, each of which processes one read (batch size 1) at a
	// time.  Each read's tensor maths should therefore stay SINGLE-THREADED,
	// otherwise LibTorch's intra-op pool (which defaults to the core count)
	// multiplies with the OpenMP pool and massively oversubscribes the CPU
	// (threads x intra_op threads).  For a high-volume workload, per-read
	// single-threaded ops + read-level parallelism gives the best throughput.
	//
	// Power users can override the intra-op count with DNASCENT_TORCH_THREADS
	// (e.g. when running with a small -t on long CPU-only reads).
	int intra_op = 1;
	const char *env = std::getenv("DNASCENT_TORCH_THREADS");
	if (env != nullptr){
		int v = std::atoi(env);
		if (v >= 1) intra_op = v;
	}
	(void) threads;   // read-level parallelism is handled by OpenMP in detect_main
	try{
		at::set_num_threads(intra_op);
	}
	catch(const std::exception &e){
		// set_num_threads throws if called after the pool is initialised;
		// safe to ignore on subsequent loads.
	}

	std::shared_ptr<ModelSession> ms = std::make_shared<ModelSession>();
	ms -> device = device;

	try{
		ms -> module = torch::jit::load(model_path, device);
	}
	catch(const c10::Error &e){
		std::cerr << "Error: failed to load TorchScript model from '" << model_path << "'." << std::endl;
		std::cerr << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}

	// Ensure eval mode (disables dropout / selects the sequential scan path).
	ms -> module.eval();

	// DNAscent feeds the model ONE read at a time and every read has a different
	// length.  The default TorchScript *profiling* graph executor re-specialises
	// on each new input shape and caches a fresh optimised graph per shape, which
	// leaks memory and makes each successive read progressively slower.  Switch
	// to the legacy executor, which specialises only on dtype/rank (not concrete
	// sizes), so every read — whatever its length — reuses the same graph.
	torch::jit::getProfilingMode()  = false;
	torch::jit::getExecutorMode()   = false;

	return ms;
}


std::shared_ptr<ModelSession> model_load_cpu(const char *model_path, unsigned int threads){

	return load_module(model_path, torch::Device(torch::kCPU), threads);
}


std::shared_ptr<ModelSession> model_load_gpu(const char *model_path, unsigned char device, unsigned int threads){

	if (not torch::cuda::is_available()){
		std::cerr << "Error: GPU inference requested but no CUDA device is available." << std::endl;
		exit(EXIT_FAILURE);
	}

	torch::Device cudaDevice(torch::kCUDA, static_cast<c10::DeviceIndex>(device));
	return load_module(model_path, cudaDevice, threads);
}
