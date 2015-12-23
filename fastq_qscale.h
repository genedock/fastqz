#include <stdbool.h>
#include "IO_stream.h"
/* The scales here are taken from the wikipedia article for FASTQ, the relavent
 * part is reproduced here:
 *
 *  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
 *  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
 *  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
 *  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
 *  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
 *  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 *  |                         |    |        |                              |                     |
 * 33                        59   64       73                            104                   126
 *
 *    S - Sanger        Phred+33,  raw reads typically (0, 40)
 *    X - Solexa        Solexa+64, raw reads typically (-5, 40)
 *    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 *    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
 *       with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
 *       (Note: See discussion above).
 *        L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 */


typedef struct qual_scale_t_
{
    const char* description;
    char min_qual, max_qual;
} qual_scale_t;

/* When the scale is ambiguous, we choose the first compatible one. Hence these
 * are ordered roughly by increasing exclusivity. */
#define NUM_SCALES 5
static qual_scale_t scales[NUM_SCALES] =
{
    {"Sanger/Phred+33",       '!', 'I'},
    {"Illumina 1.8/Phred+33", '!', 'J'},
    {"Illumina 1.5/Phred+64", 'B', 'h'},
    {"Illumina 1.3/Phred+64", '@', 'h'},
    {"Solexa/Solexa+64",      ';', 'h'}
};

static inline int fastq_qualscale(const char* fn, gzFile fin);
static inline char *readNextNode(gzFile fq,char *buf);
static inline uint32_t make_bitset(char min_qual, char max_qual);
static inline bool single_bit(uint32_t x);

/* Return true if x has excatly one 1 bit. */
bool single_bit(uint32_t x)
{
    return x && !(x & (x - 1));
}

/* Make a bitset of compatible scales. */
uint32_t make_bitset(char min_qual, char max_qual)
{
    uint32_t s = 0;
    uint32_t i;
    for (i = 0; i < NUM_SCALES; ++i) {
        if (scales[i].min_qual <= min_qual &&
            scales[i].max_qual >= max_qual) {
            s |= 1 << i;
        }
    }

    return s;
}

char *readNextNode(gzFile fq,char *buf) {
    buf=gzgets(fq,buf,1024*sizeof(char));
    if (!gzeof(fq)) {
        buf=gzgets(fq,buf,1024*sizeof(char));
        buf=gzgets(fq,buf,1024*sizeof(char));
        buf=gzgets(fq,buf,1024*sizeof(char));
        *(buf+strlen(buf)-1)='\0';
        char *quality = strdup(buf);
        return quality;
    } else {
        return NULL;
    }
}

int fastq_qualscale(const char* fn, gzFile fin)
{
    char min_qual = '~', max_qual = '!';
    char *buf=(char *)calloc(1024,sizeof(char));


    char* quality = NULL;

    /* Scales compatible with the data so far. */
    uint32_t compat_scales;

    size_t n = 100000;

    while (n-- && (quality = readNextNode(fin, buf)) !=NULL && !single_bit(compat_scales)) {
        char* q = quality;
        while (*q) {
            if (*q < min_qual) min_qual = *q;
            if (*q > max_qual) max_qual = *q;
            ++q;
        }
        free(quality);
        quality = NULL;
        compat_scales = make_bitset(min_qual, max_qual);
        if (compat_scales == 0 || single_bit(compat_scales)) break;
    }
    free(buf);
    gzrewind(fin);

    if (compat_scales == 0) {
        printf("%s: Unknown scale ['%c', '%c']\n", fn, min_qual, max_qual);
        return -1;
    }
    else {
        /* low order bit */
        unsigned int i;
        for (i = 0; !(compat_scales & (1 << i)); ++i) {}
        printf("%s: %s\n", fn, scales[i].description);
        return i;
    }
}