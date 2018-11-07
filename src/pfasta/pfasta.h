#ifndef PFASTA_H
#define PFASTA_H

#ifdef __cplusplus
extern "C" {
#endif

/** The following is the maximum length of an error string. It has to be
 * carefully chosen, so that all calls to PF_FAIL_STR succeed. For instance,
 * the line number can account for up to 20 characters.
 */
#define PF_ERROR_STRING_LENGTH 100

typedef struct pfasta_file {
	char *buffer, *readptr, *fillptr;
	char *errstr;
	int errno__;
	int fd;
	size_t line;
	char errstr_buf[PF_ERROR_STRING_LENGTH];
	char unexpected_char;
} pfasta_file;

typedef struct pfasta_seq {
	char *name, *comment, *seq;
} pfasta_seq;

int pfasta_parse(pfasta_file *, int file_descriptor);
void pfasta_free(pfasta_file *);
void pfasta_seq_free(pfasta_seq *);
int pfasta_read(pfasta_file *, pfasta_seq *);

const char *pfasta_strerror(const pfasta_file *);

#ifdef __cplusplus
}
#endif

#endif /* PFASTA_H */
