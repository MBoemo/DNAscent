#include "../src/tensor_config.h"
#include <tensorflow/c/c_api.h>

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

static void check(bool condition, const char* message){
	if (!condition){
		std::cerr << message << '\n';
		std::exit(EXIT_FAILURE);
	}
}

int main(){
	// Golden messages check values and nested lengths independently of the encoder.
	check(tensor_config::configProto(128, 128, 2) == std::vector<uint8_t>({
		0x0a, 0x08, 0x0a, 0x03, 'C', 'P', 'U', 0x10, 0x80, 0x01,
		0x10, 0x80, 0x01, 0x28, 0x02, 0x38, 0x01}), "CPU 256-thread configuration differs");
	check(tensor_config::configProto(256, 1, 256, '0') == std::vector<uint8_t>({
		0x0a, 0x08, 0x0a, 0x03, 'C', 'P', 'U', 0x10, 0x80, 0x02,
		0x0a, 0x07, 0x0a, 0x03, 'G', 'P', 'U', 0x10, 0x01,
		0x10, 0x01, 0x28, 0x80, 0x02, 0x32, 0x05, 0x20, 0x01,
		0x2a, 0x01, '0', 0x38, 0x01}), "GPU 256-thread configuration differs");
	check(tensor_config::configProto(16384, 16384, 2) == std::vector<uint8_t>({
		0x0a, 0x09, 0x0a, 0x03, 'C', 'P', 'U', 0x10, 0x80, 0x80, 0x01,
		0x10, 0x80, 0x80, 0x01, 0x28, 0x02, 0x38, 0x01}), "Three-byte varint differs");

	TF_Status* status = TF_NewStatus();
	const unsigned int threadCounts[] = {
		0, 1, 4, 63, 64, 127, 128, 200, 255, 256, 300, 511, 512,
		16383, 16384, 32767, 32768, static_cast<unsigned int>(std::numeric_limits<int>::max())
	};
	const unsigned char devices[] = {0, '0', '1', '9'};
	for (unsigned int threads : threadCounts){
		for (unsigned char device : devices){
			const unsigned int half = std::max(1u, threads/2);
			const auto config = device == 0 ? tensor_config::configProto(half, half, 2)
			                               : tensor_config::configProto(threads, 1, threads, device);
			TF_SessionOptions* options = TF_NewSessionOptions();
			TF_SetConfig(options, config.data(), config.size(), status);
			check(TF_GetCode(status) == TF_OK, TF_Message(status));
			TF_DeleteSessionOptions(options);

			// Keep the existing byte representation for all previously valid values.
			if ((device == 0 && half < 128) || (device != 0 && threads < 128)){
				const uint8_t cpus = device == 0 ? half : threads;
				const uint8_t intra = device == 0 ? half : 1;
				const uint8_t inter = device == 0 ? 2 : threads;
				const std::vector<uint8_t> legacy = device == 0
					? std::vector<uint8_t>{0x0a, 0x07, 0x0a, 0x03, 'C', 'P', 'U', 0x10, cpus, 0x10, intra, 0x28, inter, 0x38, 0x01}
					: std::vector<uint8_t>{0x0a, 0x07, 0x0a, 0x03, 'C', 'P', 'U', 0x10, cpus, 0x0a, 0x07, 0x0a, 0x03, 'G', 'P', 'U', 0x10, 0x01, 0x10, intra, 0x28, inter, 0x32, 0x05, 0x20, 0x01, 0x2a, 0x01, device, 0x38, 0x01};
				check(config == legacy, "Previously valid configuration changed");
			}

			// Machine-readable output also permits verification with ConfigProto.
			std::cout << std::dec << threads << ' ' << static_cast<unsigned int>(device) << ' ';
			for (uint8_t byte : config){
				std::cout << std::hex << std::setfill('0') << std::setw(2) << static_cast<unsigned int>(byte);
			}
			std::cout << '\n';
		}
	}
	TF_DeleteStatus(status);
}
