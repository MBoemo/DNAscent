//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include "tensor.h"
#include <algorithm>


std::pair< std::shared_ptr<ModelSession>, std::shared_ptr<TF_Graph *> > model_load_cpu_twoInputs(const char *saved_model_dir, unsigned int threads){

	int vis = setenv("CUDA_VISIBLE_DEVICES", "", 1);
	if (vis == -1){
		std::cerr << "Suppression of GPU devices failed." << std::endl;
	}

	//for CPU-only useage, the tensorflow gpu library will still print out warnings about not finding GPU/CUDA - suppress them here
	int env = setenv("TF_CPP_MIN_LOG_LEVEL", "3", 1);
	if (env == -1){
		std::cerr << "Suppression of Tensorflow logs and warnings failed." << std::endl;
	}

	std::shared_ptr<ModelSession> ms = std::make_shared<ModelSession>();
	std::shared_ptr<TF_Graph *> Graph = std::make_shared<TF_Graph *>(TF_NewGraph());

	TF_Status* Status = TF_NewStatus();

	TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
	TF_Buffer* RunOpts = NULL;

	//configure the model
	CStatus status;
	uint8_t intra_op_parallelism_threads = std::max((unsigned int)1,threads/2);
	uint8_t inter_op_parallelism_threads = 2;
	uint8_t cpus = std::max((unsigned int)1,threads/2);
	uint8_t buf[]={0xa, 0x7, 0xa, 0x3, 0x43, 0x50, 0x55, 0x10, cpus, 0x10, intra_op_parallelism_threads, 0x28, inter_op_parallelism_threads, 0x38, 0x1};
	TF_SetConfig(SessionOpts, buf,sizeof(buf),status.ptr);
	
	if(status.failure()){
		std::cout << "Model configuration failed." << std::endl;
	}

	const char* tags = "serve";
	int ntags = 1;

	std::shared_ptr<TF_Session*> session = std::make_shared<TF_Session*>(TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, saved_model_dir, &tags, ntags, *(Graph.get()), NULL, Status));
	if(TF_GetCode(Status) != TF_OK){
		printf("%s",TF_Message(Status));
	}

	auto output_op = TF_GraphOperationByName(*(Graph.get()), "StatefulPartitionedCall");
	if(!output_op){
		std::cout << "bad output name" << std::endl;
		exit(0);
	}

	ms -> session = session;
	ms -> outputs = {output_op, 0};

	return std::make_pair(ms,Graph);
}


std::pair< std::shared_ptr<ModelSession>, std::shared_ptr<TF_Graph *> > model_load_gpu_twoInputs(const char *saved_model_dir, unsigned char device, unsigned int threads){

	std::shared_ptr<ModelSession> ms = std::make_shared<ModelSession>();
	std::shared_ptr<TF_Graph *> Graph = std::make_shared<TF_Graph *>(TF_NewGraph());

	TF_Status* Status = TF_NewStatus();

	TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
	TF_Buffer* RunOpts = NULL;

	//configure the model
	CStatus status;
	uint8_t intra_op_parallelism_threads = 1;
	uint8_t inter_op_parallelism_threads = threads;
	uint8_t cpus = threads;
	uint8_t buf[]={0xa, 0x7, 0xa, 0x3, 0x43, 0x50, 0x55, 0x10, cpus, 0xa, 0x7, 0xa, 0x3, 0x47, 0x50, 0x55, 0x10, 0x1, 0x10, intra_op_parallelism_threads, 0x28, inter_op_parallelism_threads, 0x32, 0x5, 0x20, 0x1, 0x2a, 0x1, (unsigned char) device, 0x38, 0x1};
	TF_SetConfig(SessionOpts, buf,sizeof(buf),status.ptr);
	
	if(status.failure()){
		std::cout << "Model configuration failed." << std::endl;
	}

	const char* tags = "serve";
	int ntags = 1;

	std::shared_ptr<TF_Session*> session = std::make_shared<TF_Session*>(TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, saved_model_dir, &tags, ntags, *(Graph.get()), NULL, Status));
	if(TF_GetCode(Status) != TF_OK){
		printf("%s",TF_Message(Status));
	}

	auto output_op = TF_GraphOperationByName(*(Graph.get()), "StatefulPartitionedCall");
	if(!output_op){
		std::cout << "bad output name" << std::endl;
		exit(0);
	}

	ms -> session = session;
	ms -> outputs = {output_op, 0};

	return std::make_pair(ms,Graph);
}


std::shared_ptr<ModelSession> model_load_cpu(const char *saved_model_dir, unsigned int threads, const char *input_layer_name){

	std::shared_ptr<ModelSession> ms = std::make_shared<ModelSession>();
	std::shared_ptr<TF_Graph *> Graph = std::make_shared<TF_Graph *>(TF_NewGraph());

	TF_Status* Status = TF_NewStatus();

	TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
	TF_Buffer* RunOpts = NULL;

	//configure the model
	CStatus status;
	uint8_t intra_op_parallelism_threads = std::max((unsigned int)1,threads/2);
	uint8_t inter_op_parallelism_threads = 2;
	uint8_t cpus = std::max((unsigned int)1,threads/2);
	uint8_t buf[]={0xa, 0x7, 0xa, 0x3, 0x43, 0x50, 0x55, 0x10, cpus, 0x10, intra_op_parallelism_threads, 0x28, inter_op_parallelism_threads, 0x38, 0x1};
	TF_SetConfig(SessionOpts, buf,sizeof(buf),status.ptr);
	
	if(status.failure()){
		std::cout << "Model configuration failed." << std::endl;
	}

	const char* tags = "serve"; // default model serving tag; can change in future
	int ntags = 1;

	std::shared_ptr<TF_Session*> session = std::make_shared<TF_Session*>(TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, saved_model_dir, &tags, ntags, *(Graph.get()), NULL, Status));
	if(TF_GetCode(Status) != TF_OK){
		printf("%s",TF_Message(Status));
	}

	auto input_op = TF_GraphOperationByName(*(Graph.get()), input_layer_name);
	auto output_op = TF_GraphOperationByName(*(Graph.get()), "StatefulPartitionedCall");
	if(!output_op){
		std::cout << "bad output name" << std::endl;
		return nullptr;
	}
	if(!input_op){
		std::cout << "bad input name" << std::endl;
		return nullptr;
	}

	ms -> session = session;

	ms -> inputs = {input_op, 0};
	ms -> outputs = {output_op, 0};

	return ms;
}


std::shared_ptr<ModelSession> model_load_gpu(const char *saved_model_dir, unsigned char device, unsigned int threads,const char *input_layer_name){

	std::shared_ptr<ModelSession> ms = std::make_shared<ModelSession>();
	std::shared_ptr<TF_Graph *> Graph = std::make_shared<TF_Graph *>(TF_NewGraph());

	TF_Status* Status = TF_NewStatus();

	TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
	TF_Buffer* RunOpts = NULL;

	//configure the model
	CStatus status;
	uint8_t intra_op_parallelism_threads = 1;
	uint8_t inter_op_parallelism_threads = threads;
	uint8_t cpus = threads;
	uint8_t buf[]={0xa, 0x7, 0xa, 0x3, 0x43, 0x50, 0x55, 0x10, cpus, 0xa, 0x7, 0xa, 0x3, 0x47, 0x50, 0x55, 0x10, 0x1, 0x10, intra_op_parallelism_threads, 0x28, inter_op_parallelism_threads, 0x32, 0x5, 0x20, 0x1, 0x2a, 0x1, (unsigned char) device, 0x38, 0x1};
	TF_SetConfig(SessionOpts, buf,sizeof(buf),status.ptr);
	
	if(status.failure()){
		std::cout << "Model configuration failed." << std::endl;
	}

	const char* tags = "serve"; // default model serving tag; can change in future
	int ntags = 1;

	std::shared_ptr<TF_Session*> session = std::make_shared<TF_Session*>(TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, saved_model_dir, &tags, ntags, *(Graph.get()), NULL, Status));
	if(TF_GetCode(Status) != TF_OK){
		printf("%s",TF_Message(Status));
	}

	auto input_op = TF_GraphOperationByName(*(Graph.get()), input_layer_name);
	auto output_op = TF_GraphOperationByName(*(Graph.get()), "StatefulPartitionedCall");
	if(!output_op){
		std::cout << "bad output name" << std::endl;
		return nullptr;
	}
	if(!input_op){
		std::cout << "bad input name" << std::endl;
		return nullptr;
	}

	ms -> session = session;

	ms -> inputs = {input_op, 0};
	ms -> outputs = {output_op, 0};

	return ms;
}

