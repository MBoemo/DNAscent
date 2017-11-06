#pragma once
#ifndef SCRAPPIE_COMMON_H
#define SCRAPPIE_COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <assert.h>
#include <string.h>

void trim_and_segment_raw( double *raw, size_t raw_size, unsigned int *startIndex, unsigned int *endIndex );

#ifdef __cplusplus
}
#endif

#endif /* SCRAPPIE_COMMON_H */
