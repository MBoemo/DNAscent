// Too many calls to quantile lead to high failure rate for `scrappie raw`
#define BANANA 1
#include <assert.h>
#include <err.h>
#include <math.h>
#include "util.h"

int floatcmp(const void *x, const void *y) {
	float d = *(float *)x - *(float *)y;
	if (d > 0) {
		return 1;
	}
	return -1;
}


/**  Quantiles from n array
 *
 *  Using a relatively inefficent qsort resulting in O(n log n)
 *  performance but better performance is possible for small np.
 *  The array p is modified inplace, containing which quantiles to
 *  calculation on input and the quantiles on output; on error, p
 *  is filled with the value NAN.
 *
 *  @param x An array to calculate quantiles from
 *  @param nx Length of array x
 *  @param p An array containing quantiles to calculate [in/out]
 *  @param np Length of array p
 *
 *  @return void
 **/
void quantilef(const float *x, size_t nx, float *p, size_t np) {
	if (NULL == p) {
		return;
	}
	for (int i = 0; i < np; i++) {
		assert(p[i] >= 0.0f && p[i] <= 1.0f);
	}
	if (NULL == x) {
		for (int i = 0; i < np; i++) {
			p[i] = NAN;
		}
		return;
	}
	// Sort array
	float *space = malloc(nx * sizeof(float));
	if (NULL == space) {
		for (int i = 0; i < np; i++) {
			p[i] = NAN;
		}
		return;
	}
	memcpy(space, x, nx * sizeof(float));
	qsort(space, nx, sizeof(float), floatcmp);

	// Extract quantiles
	for (int i = 0; i < np; i++) {
		const size_t idx = p[i] * (nx - 1);
		const float remf = p[i] * (nx - 1) - idx;
		if (idx < nx - 1) {
			p[i] = (1.0 - remf) * space[idx] + remf * space[idx + 1];
		}
		else {
			// Should only occur when p is exactly 1.0
			p[i] = space[idx];
		}
	}

	free(space);
	return;
}


/** Median of an array
 *
 *  Using a relatively inefficent qsort resulting in O(n log n)
 *  performance but O(n) is possible.
 *
 *  @param x An array to calculate median of
 *  @param n Length of array
 *
 *  @return Median of array on success, NAN otherwise.
 **/
float medianf(const float *x, size_t n) {
	float p = 0.5;
	quantilef(x, n, &p, 1);
	return p;
}


/** Median Absolute Deviation of an array
 *
 *  @param x An array to calculate the MAD of
 *  @param n Length of array
 *  @param med Median of the array.  If NAN then median is calculated.
 *
 *  @return MAD of array on success, NAN otherwise.
 **/
float madf(const float *x, size_t n, const float *med) {
	const float mad_scaling_factor = 1.4826;
	if (NULL == x) {
		return NAN;
	}
	if (1 == n) {
		return 0.0f;
	}

	float *absdiff = malloc(n * sizeof(float));
	if (NULL == absdiff) {
		return NAN;
	}

	const float _med = (NULL == med) ? medianf(x, n) : *med;

	for (size_t i = 0; i < n; i++) {
		absdiff[i] = fabsf(x[i] - _med);
	}

	const float mad = medianf(absdiff, n);
	free(absdiff);
	return mad * mad_scaling_factor;
}
