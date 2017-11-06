#include "scrappie_common.h"
#include "util.h"


/**  Simple segmentation of a raw read by thresholding the MAD
 *
 *  The MAD of the raw signal is calculated for non-overlapping chunks and then
 *  thresholded to find regions at the beginning and end of the signal that have
 *  unusually low variation (generally a stall or open pore).  The threshhold is
 *  derived from the distribution of the calaculated MADs.
 *
 *  The threshold is chosen to be high since a single chunk above it will trigger
 *  the end of the trimming: the threshhold is chosen so it is unlikely to be
 *  exceeded in the leader but commonly exceeded in the main read.
 *
 *  @param rt Structure containing raw signal
 *  @param varseg_chunk Size of non-overlapping chunks
 *  @param varseg_thresh  The quantile to be calculated to use for threshholding
 *
 *  @return A range structure containing new start and end for read
 **/
void trim_and_segment_raw( double *raw, size_t raw_size, unsigned int *startIndex, unsigned int *endIndex ) {

	int trim_start = 200;
	int trim_end = 10;
	int varseg_chunk = 100;
	float varseg_thresh = 0.0f;

	unsigned int rawStart = 0;
	unsigned int rawEnd = raw_size;

	assert(varseg_chunk > 1);
	assert(varseg_thresh >= 0.0 && varseg_thresh <= 1.0);

	const size_t nchunk = raw_size / varseg_chunk;
	// Truncation of end to be consistent with Sloika
	rawEnd = nchunk * varseg_chunk;

	float *madarr = malloc(nchunk * sizeof(float));

	for (size_t i = 0; i < nchunk; i++) {
		madarr[i] = madf(raw + rawStart + i * varseg_chunk, varseg_chunk, NULL);
	}
	quantilef(madarr, nchunk, &varseg_thresh, 1);

	const float thresh = varseg_thresh;
	for (size_t i = 0; i < nchunk; i++) {
		if (madarr[i] > thresh) {
			break;
		}
		rawStart += varseg_chunk;
	}
	for (size_t i = nchunk; i > 0; i--) {
		if (madarr[i - 1] > thresh) {
			break;
		}
		rawEnd -= varseg_chunk;
	}
	assert(rawEnd > rawStart);

	free(madarr);

	rawStart += trim_start;
	rawEnd -= trim_end;

	*startIndex = rawStart;
	*endIndex = rawEnd;

	if (rawStart >= rawEnd) {
		free( raw );
	}
}
