#ifndef EVENT_DETECTION_H
#define EVENT_DETECTION_H

#include <stddef.h>

typedef struct {
    size_t window_length1;
    size_t window_length2;
    float threshold1;
    float threshold2;
    float peak_height;
} detector_param;


static detector_param const event_detection_defaults = {
    .window_length1 = 3,
    .window_length2 = 6,
    .threshold1 = 1.4f,
    .threshold2 = 9.0f,
    .peak_height = 0.2f
};

typedef struct {
	int pos;
	int state;
	unsigned int start;
	float length;
	float mean;
	float stdv;
} event_s;



event_s * detect_events(int *raw, size_t raw_size, detector_param const edparam);

#endif    
