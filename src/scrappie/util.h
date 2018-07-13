#pragma once
#ifndef UTIL_H
#define UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <immintrin.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

void quantilef(const float *x, size_t nx, float *p, size_t np);
float medianf(const float *x, size_t n);
float madf(const float *x, size_t n, const float *med);

#ifdef __cplusplus
}
#endif

#endif                          /* UTIL_H */
