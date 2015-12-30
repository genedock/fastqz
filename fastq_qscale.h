#include <stdbool.h>
#include <stdint.h>
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

typedef struct _fastq_ {
    char *name;
    char *seq;
    char *quality;
} Fastq;

#define _putc(__ch, __out) *__out++ = (__ch)
#define _getc(in, in_) (in<in_?(*in++):-1)
#define _rewind(in,_in) in = _in

#define free_Fastq(line)    \
do{                         \
    free((line)->name);     \
    free((line)->seq);      \
    free((line)->quality);  \
}while (0)

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
static inline char *readNext(gzFile fq,char *buf);
static inline uint32_t make_bitset(char min_qual, char max_qual);
static inline bool single_bit(uint32_t x);
static inline int readNextNode(gzFile fq,char *buf,Fastq *line);
static inline int getNumCores();

#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

int getNumCores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

void int2str(int c, int base, char *ret)
{
    const char *tab = "0123456789abcdef";
    if (c == 0) ret[0] = '0', ret[1] = 0;
    else {
        int l=0, x=0, y=0;
        char buf[16];
        for (l = 0, x = c < 0? -c : c; x > 0; ) {
            buf[l++] = tab[x%base];
            x /= base;
        }
        if (c < 0) buf[l++] = '-';
        for (x = l - 1, y = 0; x >= 0; --x) ret[y++] = buf[x];
        ret[y] = 0;
    }
}

/* Return true if x has excatly one 1 bit. */
bool single_bit(uint32_t x) {
    return x && !(x & (x - 1));
}

/* Make a bitset of compatible scales. */
uint32_t make_bitset(char min_qual, char max_qual) {
    uint32_t s = 0;
    uint32_t i;
    for (i = 0; i < NUM_SCALES; ++i) {
        if (scales[i].min_qual <= min_qual && scales[i].max_qual >= max_qual) {
            s |= 1 << i;
        }
    }
    return s;
}

char *readNext(gzFile fq,char *buf) {
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

char *readNextFile(FILE *fq,char *buf) {
    buf=fgets(buf,1024*sizeof(char),fq);
    if (!feof(fq)) {
        buf=fgets(buf,1024*sizeof(char),fq);
        buf=fgets(buf,1024*sizeof(char),fq);
        buf=fgets(buf,1024*sizeof(char),fq);
        *(buf+strlen(buf)-1)='\0';
        char *quality = strdup(buf);
        return quality;
    } else {
        return NULL;
    }
}

int readNextNode(gzFile fq,char *buf,Fastq *line) {
    buf=gzgets(fq,buf,1024*sizeof(char));
    if (!gzeof(fq)) {
        *(buf+strlen(buf)-1)=0; //remove the last \n in the buf
        line->name  = strdup(buf);
        
        buf=gzgets(fq,buf,1024*sizeof(char));
        *(buf+strlen(buf)-1)=0; //remove the last \n in the buf
        line->seq = strdup(buf);

        buf=gzgets(fq,buf,1024*sizeof(char));

        buf=gzgets(fq,buf,1024*sizeof(char));
        line->quality = strdup(buf);
        return 0;
    } else {
        return 1;
    }
}

int readNextNodeFile(FILE *fq,char *buf,Fastq *line) {
    buf=fgets(buf,1024*sizeof(char),fq);
    if (!feof(fq)) {
        *(buf+strlen(buf)-1)=0; //remove the last \n in the buf
        line->name  = strdup(buf);
        
        buf=fgets(buf,1024*sizeof(char),fq);
        *(buf+strlen(buf)-1)=0; //remove the last \n in the buf
        line->seq = strdup(buf);

        buf=fgets(buf,1024*sizeof(char),fq);

        buf=fgets(buf,1024*sizeof(char),fq);
        line->quality = strdup(buf);
        return 0;
    } else {
        return 1;
    }
}

int fastq_qualscaleFile(const char* fn, FILE *fin) {
    char min_qual = '~', max_qual = '!';
    char *buf=(char *)calloc(1024,sizeof(char));
    char* quality = NULL;
    /* Scales compatible with the data so far. */
    uint32_t compat_scales=0;
    size_t n = 100000;
    while (n-- && (quality = readNextFile(fin, buf)) !=NULL && !single_bit(compat_scales)) {
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
    rewind(fin);

    if (compat_scales == 0) {
        fprintf(stderr,"%s: Unknown scale ['%c', '%c']\n", fn, min_qual, max_qual);
        return -1;
    }
    else {
        /* low order bit */
        unsigned int i;
        for (i = 0; !(compat_scales & (1 << i)); ++i) {}
        fprintf(stderr,"%s: %s\n", fn, scales[i].description);
        return i;
    }
}

int fastq_qualscale(const char* fn, gzFile fin) {
    char min_qual = '~', max_qual = '!';
    char *buf=(char *)calloc(1024,sizeof(char));
    char* quality = NULL;
    /* Scales compatible with the data so far. */
    uint32_t compat_scales;
    size_t n = 100000;
    while (n-- && (quality = readNext(fin, buf)) !=NULL && !single_bit(compat_scales)) {
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
        fprintf(stderr,"%s: Unknown scale ['%c', '%c']\n", fn, min_qual, max_qual);
        return -1;
    }
    else {
        /* low order bit */
        unsigned int i;
        for (i = 0; !(compat_scales & (1 << i)); ++i) {}
        fprintf(stderr,"%s: %s\n", fn, scales[i].description);
        return i;
    }
}
