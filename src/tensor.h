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
#include <cassert>
#include <fstream>
#include "../tensorflow/include/tensorflow/c/eager/c_api.h"

#define MAX_DIM 16

//start: adapted from https://github.com/aljabr0/from-keras-to-c
//licensed under Apache 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
class CStatus{
	public:
		TF_Status *ptr;
		CStatus(){
			ptr = TF_NewStatus();
		}

		void dump_error()const{
			std::cerr << "TF status error: " << TF_Message(ptr) << std::endl;
		}

		inline bool failure()const{
			return TF_GetCode(ptr) != TF_OK;
		}

		~CStatus(){
			if(ptr)TF_DeleteStatus(ptr);
		}
};

namespace detail {

	template<class T>
	class TFObjDeallocator;

	template<>
	struct TFObjDeallocator<TF_Status> { static void run(TF_Status *obj) { TF_DeleteStatus(obj); }};

	template<>
	struct TFObjDeallocator<TF_Graph> { static void run(TF_Graph *obj) { TF_DeleteGraph(obj); }};

	template<>
	struct TFObjDeallocator<TF_Tensor> { static void run(TF_Tensor *obj) { TF_DeleteTensor(obj); }};

	template<>
	struct TFObjDeallocator<TF_SessionOptions> { static void run(TF_SessionOptions *obj) { TF_DeleteSessionOptions(obj); }};

	template<>
	struct TFObjDeallocator<TF_Buffer> { static void run(TF_Buffer *obj) { TF_DeleteBuffer(obj); }};

	template<>
	struct TFObjDeallocator<TF_ImportGraphDefOptions> {
		static void run(TF_ImportGraphDefOptions *obj) { TF_DeleteImportGraphDefOptions(obj); }
	};

	template<>
	struct TFObjDeallocator<TF_Session> {
		static void run(TF_Session *obj) {
			CStatus status;
			TF_DeleteSession(obj, status.ptr);
			if (status.failure()) {
				status.dump_error();
			}
		}
	};
}


template<class T> struct TFObjDeleter{
	void operator()(T* ptr) const{
		detail::TFObjDeallocator<T>::run(ptr);
	}
};


template<class T> struct TFObjMeta{
	typedef std::unique_ptr<T, TFObjDeleter<T>> UniquePtr;
};


template<class T> typename TFObjMeta<T>::UniquePtr tf_obj_unique_ptr(T *obj){
	typename TFObjMeta<T>::UniquePtr ptr(obj);
	return ptr;
};


class ModelSession{
	public:
		std::shared_ptr<TF_Graph*> graph;
		std::shared_ptr<TF_Session*> session;
		TF_Output inputs, outputs;
};


template<class T> static void free_cpp_array(void* data, size_t length){
	delete []((T *)data);
}


template<class T> static void cpp_array_deallocator(void* data, size_t length, void* arg){
	delete []((T *)data);
}


struct TensorShape{
	int64_t values[MAX_DIM];
	int dim;

	int64_t size()const{
		assert(dim>=0);
		int64_t v=1;
		for(int i=0; i<dim; i++) v*=values[i];
		return v;
    }
};

//end adapted from https://github.com/aljabr0/from-keras-to-c


std::shared_ptr<ModelSession> model_load_cpu(const char *filename, unsigned int threads, const char *);
std::shared_ptr<ModelSession> model_load_gpu(const char *filename, unsigned char device, unsigned int threads, const char *);
std::pair< std::shared_ptr<ModelSession>, std::shared_ptr<TF_Graph *> > model_load_cpu_twoInputs(const char *filename, unsigned int threads);
std::pair< std::shared_ptr<ModelSession>, std::shared_ptr<TF_Graph *> > model_load_gpu_twoInputs(const char *filename, unsigned char device, unsigned int threads);

#endif
