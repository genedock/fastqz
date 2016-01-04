//gcc -g -O3 -Wall gzfastq_sort_list.c  -o gzfastq_sort_list -I/share/software/software/zlib_1.2.8_install/include -I. -L. -L/share/software/software/zlib_1.2.8_install/lib -lz -llist
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <err.h>
#include <time.h>
#include <sys/time.h>
#include "fastq_qscale.h"

struct globalArgs_t {
	const char *infile;
	const char *outfile;
} globalArgs;

void display_usage(char * argv[]);
void load_file(FILE *fp,const char *outfile);

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(10240*sizeof(char));
	const char* usage=
"\nCopyright (c)  20155555\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Discription:\n  This program is used for sorting fastq files by name.\n" \
"Usage: %s [-i Infile] [-o OUTFILE] [-r reads_count_for allocating memory] [-s|-n] [-h] \n" \
"Example1:\n  zcat 13C37198_L7_I012.R1.clean.fastq.gz | %s -r 8000000 -o out -s\n" \
"\n" \
"   [-i Infile] = Infile.                                      [required]\n" \
"   [-o OUTPUT] = OUTPUT file.default is \'out\'.                 [option]\n" \
"   [-r reads_num] = reads_count_for allocating memory.         [option]\n" \
"   [-s ] sort by seqyence.                                     [option]\n" \
"   [-n ] sort by sequence name.                                [option]\n" \
"   [-h] This helpful help screen.                              [option]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	free(buffer);
	exit(1);
}

void load_file(FILE *fp,const char *outfile) {
	int threads=getNumCores(),blockSize=640000;
	Fastq **chunk=(Fastq **)calloc(threads,sizeof(Fastq *));
	int i,*chunkLen = (int *)calloc(threads,sizeof(int));
	for (i=0;i<threads;i++){
    	chunk[i]=(Fastq *)calloc(blockSize,sizeof(Fastq));
    }
	char *LineBuf = (char *)calloc(4096,sizeof(char));
	int eofFlag=0,line=0,j=0,k=0,blockCount=0;
	long long begin=usec();
	FILE *out= fcreat_outfile(outfile,".pigz.fq");
	for (line=0,k=0; 1; line++) {
      int blockIndex = line % blockSize;
      if ((!blockIndex && line) || eofFlag) {
        chunkLen[k++ % threads]=!eofFlag ? blockSize : blockIndex-1;
        if (!eofFlag||(eofFlag && blockIndex>1)) ++blockCount;
        int chunkIndex = blockCount % threads;
        if ((!chunkIndex && line) || eofFlag) {
          int chunkCount = !eofFlag ? threads : chunkIndex;
          for (j=0;j<chunkCount ;j++) {
            int bi=0;
            for(bi=0;bi<chunkLen[j];bi++){
            	fprintf(stdout,"%s\n%s\n+\n%s\n",chunk[j][bi].name,chunk[j][bi].seq,chunk[j][bi].quality);
            	free_Fastq(chunk[j]+bi);
            }
          }
        }
      }
      if (eofFlag) break;
      eofFlag = readNextNodeFile(fp,LineBuf,&chunk[blockCount % threads][blockIndex]);
    }
	fprintf(stderr,"done read file at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
	free(chunk);
	free(LineBuf);
	fclose(out);
	fprintf(stderr,"done free list at %.3f s\n",(double)(usec()-begin)/CLOCKS_PER_SEC);
}

int main(int argc, char *argv[]) {
	int opt = 0;
	globalArgs.infile="-";
	globalArgs.outfile="-";
	const char *optString = "i:o:h?";
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'i':
				globalArgs.infile = optarg;
				break;
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case '?':	/* fall-through is intentional */
			case 'h':
				display_usage(argv);
				break;
			default:
				fprintf(stderr,"error parameter!\n");
				break;
		}
		opt = getopt( argc, argv, optString );
	}

	FILE *fp=fopen_input_stream(globalArgs.infile);
	load_file(fp,globalArgs.outfile);
	fclose(fp);
	return 0;
}
