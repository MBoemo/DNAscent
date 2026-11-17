//----------------------------------------------------------
// This software is licensed under GPL-3.0.
//----------------------------------------------------------

#ifndef TENSOR_CONFIG_H_
#define TENSOR_CONFIG_H_

#include <cstdint>
#include <vector>

namespace tensor_config {

inline void appendVarint(std::vector<uint8_t>& buffer, uint64_t value){
	while (value >= 0x80){
		buffer.push_back(static_cast<uint8_t>((value & 0x7f) | 0x80));
		value >>= 7;
	}
	buffer.push_back(static_cast<uint8_t>(value));
}

inline void appendMessage(std::vector<uint8_t>& buffer, uint8_t tag, const std::vector<uint8_t>& message){
	buffer.push_back(tag);
	appendVarint(buffer, message.size());
	buffer.insert(buffer.end(), message.begin(), message.end());
}

// Serialize tensorflow.ConfigProto without depending on the C++ protobuf API.
// device is the existing ASCII GPU selector ('0', '1', ...); zero means CPU.
inline std::vector<uint8_t> configProto(unsigned int cpus, unsigned int intra, unsigned int inter, unsigned char device = 0){
	std::vector<uint8_t> config;
	std::vector<uint8_t> cpuCount = {0x0a, 0x03, 'C', 'P', 'U', 0x10};
	appendVarint(cpuCount, cpus);
	appendMessage(config, 0x0a, cpuCount); // device_count["CPU"]
	if (device != 0){
		appendMessage(config, 0x0a, {0x0a, 0x03, 'G', 'P', 'U', 0x10, 0x01});
	}
	config.push_back(0x10); // intra_op_parallelism_threads
	appendVarint(config, intra);
	config.push_back(0x28); // inter_op_parallelism_threads
	appendVarint(config, inter);
	if (device != 0){
		// GPUOptions: allow_growth and the one-character visible_device_list.
		appendMessage(config, 0x32, {0x20, 0x01, 0x2a, 0x01, device});
	}
	config.insert(config.end(), {0x38, 0x01}); // allow_soft_placement
	return config;
}

}

#endif
