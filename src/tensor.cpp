//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include "tensor.h"
//#include "../tensorflow/include/tensorflow/c/c_api.h"
#include <algorithm>


//start: adapted from https://github.com/aljabr0/from-keras-to-c
//licensed under Apache 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
static TF_Buffer* read_tf_buffer_from_file(const char* file) {
	std::ifstream t(file, std::ifstream::binary);
	t.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	t.seekg(0, std::ios::end);
	size_t size = t.tellg();
	auto data = std::make_unique<char[]>(size);
	t.seekg(0);
	t.read(data.get(), size);

	TF_Buffer *buf = TF_NewBuffer();
	buf->data = data.release();
	buf->length = size;
	buf->data_deallocator = free_cpp_array<char>;
	return buf;
}
//end adapted from https://github.com/aljabr0/from-keras-to-c


std::shared_ptr<ModelSession> model_load_cpu(const char *filename, const char *input_name, const char *output_name, unsigned int threads){

	int vis = setenv("CUDA_VISIBLE_DEVICES", "", 1);
	if (vis == -1){
		std::cerr << "Suppression of GPU devices failed." << std::endl;
	}

	//for CPU-only useage, the tensorflow gpu library will still print out warnings about not finding GPU/CUDA - suppress them here
	int env = setenv("TF_CPP_MIN_LOG_LEVEL", "3", 1);
	if (env == -1){
		std::cerr << "Suppression of Tensorflow logs and warnings failed." << std::endl;
	}

	CStatus status;
	std::shared_ptr<ModelSession> ms = std::make_shared<ModelSession>();
	std::shared_ptr<TF_Graph *> graph = std::make_shared<TF_Graph *>(TF_NewGraph());

	{
        // Load a protobuf containing a GraphDef
        auto graph_def=read_tf_buffer_from_file(filename);
        if(!graph_def) return nullptr;
        auto graph_opts=TF_NewImportGraphDefOptions();
        TF_GraphImportGraphDef(*(graph.get()), graph_def, graph_opts, status.ptr);
	}

	if(status.failure()){
		status.dump_error();
		return nullptr;
	}
	ms -> graph = graph;

	auto input_op = TF_GraphOperationByName(*(graph.get()), input_name);
	auto output_op = TF_GraphOperationByName(*(graph.get()), output_name);
	if(!input_op || !output_op){
		return nullptr;
	}

	TF_SessionOptions *opts = TF_NewSessionOptions();

	//set multithreading
	//the following buffer is equivalent to
	//config = tf.ConfigProto(allow_soft_placement=True,device_count = {'CPU':<threads>/2},intra_op_parallelism_threads=<threads>/2,inter_op_parallelism_threads=2)
	uint8_t intra_op_parallelism_threads = std::max((unsigned int)1,threads/2);
	uint8_t inter_op_parallelism_threads = 2;
	uint8_t cpus = std::max((unsigned int)1,threads/2);
	uint8_t buf[]={0xa, 0x7, 0xa, 0x3, 0x43, 0x50, 0x55, 0x10, cpus, 0x10, intra_op_parallelism_threads, 0x28, inter_op_parallelism_threads, 0x38, 0x1};

 	TF_SetConfig(opts, buf,sizeof(buf),status.ptr);

	std::shared_ptr<TF_Session*> session = std::make_shared<TF_Session*>(TF_NewSession(*(graph.get()), opts, status.ptr));

	if(status.failure()){
		return nullptr;
	}
	assert(session);
	ms -> session = session;

	ms -> inputs = {input_op, 0};
	ms -> outputs = {output_op, 0};

	return ms;
}


std::shared_ptr<ModelSession> model_load_gpu(const char *filename, const char *input_name, const char *output_name, unsigned char device, unsigned int threads){

	CStatus status;
	std::shared_ptr<ModelSession> ms = std::make_shared<ModelSession>();
	std::shared_ptr<TF_Graph *> graph = std::make_shared<TF_Graph *>(TF_NewGraph());

	{
        // Load a protobuf containing a GraphDef
        auto graph_def=read_tf_buffer_from_file(filename);
        if(!graph_def) return nullptr;
        auto graph_opts=TF_NewImportGraphDefOptions();
        TF_GraphImportGraphDef(*(graph.get()), graph_def, graph_opts, status.ptr);
	}

	if(status.failure()){
		status.dump_error();
		return nullptr;
	}
	ms -> graph = graph;

	auto input_op = TF_GraphOperationByName(*(graph.get()), input_name);
	auto output_op = TF_GraphOperationByName(*(graph.get()), output_name);
	if(!input_op || !output_op){
		return nullptr;
	}

	//the buffer that follows is equivalent to:
	//config = tf.ConfigProto(allow_soft_placement=True,log_device_placement=False,device_count = {'GPU': 1,'CPU':<threads>},intra_op_parallelism_threads=<threads>/2,inter_op_parallelism_threads=2)
	//config.gpu_options.allow_growth=True
	//config.gpu_options.visible_device_list= <device>
	uint8_t intra_op_parallelism_threads = 1;
	uint8_t inter_op_parallelism_threads = threads;
	uint8_t cpus = threads;
	uint8_t buf[]={0xa, 0x7, 0xa, 0x3, 0x43, 0x50, 0x55, 0x10, cpus, 0xa, 0x7, 0xa, 0x3, 0x47, 0x50, 0x55, 0x10, 0x1, 0x10, intra_op_parallelism_threads, 0x28, inter_op_parallelism_threads, 0x32, 0x5, 0x20, 0x1, 0x2a, 0x1, device, 0x38, 0x1};

	TF_SessionOptions *opts = TF_NewSessionOptions();
	TF_SetConfig(opts, buf,sizeof(buf),status.ptr);
	//TF_EnableXLACompilation(opts,true);
	std::shared_ptr<TF_Session*> session = std::make_shared<TF_Session*>(TF_NewSession(*(graph.get()), opts, status.ptr));

	if(status.failure()){
		return nullptr;
	}
	assert(session);
	ms -> session = session;

	ms -> inputs = {input_op, 0};
	ms -> outputs = {output_op, 0};

	return ms;
}

