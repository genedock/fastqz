/* fastqz v1.5 - Sanger FASTQ compressor

  Copyright (C) 2012, Matt Mahoney, Dell Inc.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

TO COMPILE

g++ -O3 -msse2 -s -lpthread fastqz.cpp libzpaq.cpp -o fastqz

You need libzpaq.cpp and libzpaq.h from either
https://sourceforge.net/projects/zpaq/ or
http://mattmahoney.net/zpaq/
libzpaq is public domain.

Also, to use in Windows you need to install Pthreads-Win32 from
http://sourceware.org/pthreads-win32/
In particular you need pthread.h to compile and pthreadGC2.dll
in your PATH to run. Pthreads-Win32 is licensed under LGPL.

libzpaq uses Just-In-Time (JIT) optimization of ZPAQL code on
x86 32 or 64 bit processors. To run on other processors, compile
with -DNOJIT to disable this feature. It will still work but run slower.


USAGE

fx is a compressor for Sanger FASTQ files. It has two compression modes,
fast and slow.

Usage: fastqz command input output [reference]
Commands:
  c - compress input to output.fx?.zpaq (3 files, ? = {h,b,q})
  d - decompress input.fx?.zpaq to output
  e - encode input to output.fx? without zpaq compression (faster)
  f - decode input.fx? to output
  cQ, eQ - quantize quality values to 35 plus a multiple of Q, rounding
           down. Default is c1 or e1.

Commands c and e compress. c compresses smaller but e compresses faster.
The corresponding decompression commands are d and f respectively.
You need 1.5 GB memory to compress with c or decompress with d.
They also both produce temporary files taking the same space as the
output of command e. The e and f commands don't use significant
memory and don't produce temporary files.

Using a quantization like c2 or e4 is lossy but improves compression
when exact quality values are not needed. Values are rounded down.

Compression produces 3 files. Command e produces files named
output.fxh, output.fxb, output.fxq. Command c produces files named
output.fxh.zpaq, output.fxb.zpaq, output.fxq.zpaq
When decompressing, omit the .fx? or .fx?.zpaq extension
on the input file names. The extensions will be assumed.

Input for compression is expected to be a Sanger FASTQ file.
The file consists of "reads" from a DNA sequencing machine. Each
read has the following format:

  @header
  ACGTN    (base calls, length n)
  +
  !..I#    (quality scores, length n, ASCII 34..73 for A,C,G,T, 33 for N)

Maximum line length is 4095. Lines must be terminated by LF
(ASCII 10) only (no CR). All base and quality lines must have
the same length (read length = n) throughout the file.
Files not in this format are rejected.

If [reference] is present, then it is the file name of a reference
genome. The same reference must be present for decompression.
The reference genome consists of a sequence of bases packed 4
per byte in MSB to LSB order with ACGT=0..3. You can use the
program fapack to convert FASTA files into this format.
The reference genome cannot be bigger than 1 GB (2^32 bases).
You need 1.5 GB memory to encode and 1 GB to decode.
A fourth file will be produced: output.fxa.zpaq or output.fxa
containing compressed alignments.


COMPRESSION FORMAT

Command "c" and "e" both split the input into 3 or 4 parts and
compress them as described below. Command "c" further compresses
each of the 3 or 4 files using a different ZPAQ model.

Headers (.fxh) are coded in the form (j,k,len,xxx...,0) which means
go to column j-1 (first column is 0) in the previous header and
add k-1 to the decimal number ending there. If k=1, then skip
this step. Then copy the first len characters of the modified previous
header, then output xxx, and finally a linefeed (ASCII 10). Save this
output, minus the linefeed.

The first 2 bytes of the .fxh file encodes the read length, n,
MSB first (e.g. 0,100 if all base and quality lines have length 100).

Base calls (.fxb) are encoded first by deleting all N's. These can be
restored because their location is indicated by a quality score
of 33. Then the remaining bases are encoded in self terminating
base 4 with A=1, T=2, C=3, G=4 allowing 3 or 4 bases per byte.
For example, "TACT" is coded as 2*64 + 1*16 + 3*4 + 2*1 = 158.

If a reference is given, then a list of matches are stored in a .fxa
file. The format is:

  (m1+1+128*dir,m2+1,m3+1,m4+1,p3,p2,p1,p0)  to encode a match
  (0)                                        to encode no match

where p3..p0 is a 32 bit pointer (MSB first)
into the reference genome after expanding to 1 base per element
(with 0..3=ACGT) and padding the ends with 16384 zeros (or A).
'dir' is 0 for a match in the forward direction or 1 for a
match in the reverse direction starting at the same point but
exchanging A with T and C with G. m1..m4 are the locations of
the first 4 diferences between the base sequence (after deleting
N's) and the reference, in the range 0..len-1 where len is the
length of the sequence with N's deleted. Thus, the bytes are
coded in the range 1..len, with bit 7 of the first byte set if
the match is reversed. The mismatches are in ascending order.
If there are less than 4 mismatches, then the remaining bytes
are coded as len+1. Thus, only reads up to 126 can be fully
matched.

If a match is present, then only the corresponding mismatched bases,
plus any bases after m4 (except N), are written to the .fxb file.
If the first byte is 0, then there is no match and the entire
base string is written (except N).

Quality scores are decoded as follows: q=1..72 decode as q+32
(33..104). q=73..136 decode as a pair (q-73)%8+64, (q-73)/8+64,
both in the range 64..71. q=137..200 decode as the triple
(q-137)%4+68, (q-137)/4%4+68, (q-137)/16+68 in the range 68..71.
q=201..255 decodes as 71 repeated q-200 (1..55) times. q=0
decodes by setting all remaining values to 35 and terminating
the sequence. The coding takes advantage of the high frequency
of q at or just below 71 that group early in the sequence, and
of sequences that end in runs of 35.

Command "c" further compresses the output.fx? files
to output.fx?.zpaq files as defined by the ZPAQ level 2 standard
which can be found at http://mattmahoney.net/zpaq/ or
https://sourceforge.net/projects/zpaq/

ZPAQ is a configurable compression format based on the PAQ context
mixing algorithm with bit-wise prediction and arithmetic coding.
Context models are described in ZPAQL byte code, which is saved to
the compressed file and can be read by a generic ZPAQ decompressor.
Thus, a FASTQ file compressed with "fastqz c" could be decompressed
first with zpaq and then with "fastq f" as opposed to decompressing
with "fastq d".

ZPAQL byte code describes an array of components and code to compute
contexts. Each component takes a context and possibly the predictions
of earlier components and outputs a new probability that the next
bit will be a 1. The output of the last component is used to arithmetic
encode or decode the next bit. After encoding or decoding, the bit
is used to update the models to reduce their prediction errors.

Whole-byte contexts are computed on byte boundaries by code running on
a ZPAQL virtual machine. This program is executed once after modeling
each byte with that byte as input. The output is saved in an array
of 32-bit values which is available as input to the array of components.
These values are combined with the previously coded bits of the current
byte to form a complete context.

A ZPAQ model is described by a config file. In this program, the
compiled byte code is fed to the model during compression, or read
from the compressed file header during decompression. The source code
for each model is given below, followed by an explanation of the code.
The command "zpaq -mfx? l" will generate the byte code used in this
program from the sources below named "fx?.cfg" (where ? is h,b,q,a).

A config file has 3 sections:

  COMP - describes the array of modeling components.
  HCOMP - ZPAQL code to compute contexts.
  POST/PCOMP - ZPAQL code for post-processing.

Post-processing is not used, so each file ends with POST 0 END.
Modeled bits are output directly.

A ZPAQL virtual machine has 32-bit registers A,B,C,D, an array
of bytes M, an array of 32 bit unsigned integers H, a condition flag F,
and a 16 bit program counter. H is the context output to the model.
A is the input byte and accumulator for arithmetic and logical operations.
B and C are pointers into M. D points to H. *B, *C, *D refer
to the elements pointed to, modulo the array sizes. The sizes are
given by the first 2 parameters after COMP.


HEADER MODELING

(fxh.cfg model to compress headers)
comp 3 8 0 0 5 (H has size 2^3, M has size 2^8)
  0 cm 20 128  (direct 20-bit context model with max count 128*4)
  1 cm 22 128
  2 icm 18     (indirect context model with 2^(18+6) bit histories)
  3 icm 19
  4 mix 13 0 4 24 255 (13 bit context, mix 0..0+4-1, rate 24, mask 255)
hcomp
  *c=a c++ a== 0 if c=0 endif (save input in buffer M pointed to by C)
  d=0 *d=0 b=c a=c hashd (context H[0] is a hash of column number)
    a=*b hashd (combined with the byte above, saved in M)
    b-- a=*b hashd (combined with the byte to the left (order 1))
  a=*d d++ *d=a b-- a=*b hashd (context H[1] as above but order 2)
  a=*d d++ *d=a b-- a=*b hashd (context H[2] as above but order 3)
  a=*d d++ *d=a b-- a=*b hashd (context H[3] as above put order 4)
  d++ a=c a<<= 8 *d=a (context H[5] for mixer is just the column number)
  halt
post 0 end (no post-processing)

The headers are compressed using a mixture of 4 context models.
The first two are direct (CM: context -> bit prediction)
and 3 and 4 are indirect (ICM: context -> bit history -> prediction).
The context for the first model is the column number, the byte
above and the byte to the left. The next 3 add 1 to 3
more bytes to the left as context, respectively. The four
bit predictions are mixed by weighted averaging in the logistic
domain (log p/(1-p)) and the weights adapted to reduce prediction
errors. The mixer weight vector is selected by a context consisting
of the column number and the previously coded bits of the
current byte. The resulting bits are arithmetic coded.

In the code above, *C=A saves the input byte in M. C++ advances
to the next byte, which was saved from the previous line.
"A== 0 IF C=0 ENDIF" tests if the input is 0, marking the end of a
header line, and if so, resets the pointer C to the beginning of
the buffer.

The next 3 lines set the context for component 0, pointed to by D.
HASHD computes the hash *D=(*D+A+512)*773.

The next 3 lines set the contexts for components 1 through 3 by
copying the previous context hash and combining it with the next
byte back in the history buffer maintained in M and pointed to
by *B.

The last line uses the low 5 bits of the column number (in C)
as part of the 13 bit context to the mixer. The low 8 bits are
left as zeros so that during modeling the bits from the partial
byte can be added.


BASE CALL MODELING

(fxb.cfg model to compress base calls)
comp 3 3 0 0 7 (hh hm ph pm n)
  0 cm 9 255 (2 KB)
  1 cm 18 255 (1 MB)
  2 cm 25 255 (128 MB)
  3 icm 22 (256 MB)
  4 isse 23 3 (512 MB)
  5 match 26 28 (256 MB hash table, 256 MB buffer)
  6 mix 8 0 6 12 255 (order 0 mix of 0..0+6-1, rate 12, mask 255)
hcomp
  c++ *c=a b=c a=0 (save in rotating buffer M)
  d= 1 hash *d=a
  b-- d++ hash *d=a
  b-- d++ hash *d=a
  b-- d++ hash *d=a
  b-- d++ hash *d=a
  halt
post
  0
end

Base calls are modeled using an order 0..5 mix. Orders 0, 1, and 2
are direct, slow adapting (rate = error/count up to 255*4) context models.
Order 3 is indirect. Order 4 is indirect and chained to the order 3
output, i.e. order 3 prediction is mixed with a constant 1 in the
logistic domain by a pair of adaptive weights selected by the
bit history indexed by the order 4 context hash. The order 5
context is a match model which looks up the previous occurrence
of the context hash and predicts whatever bit followed. The
mixer context is bytewise order 0.

The HASH instruction computes A=(A+*B+512)*773.


QUALITY MODELING

(fxq.cfg model used to compress quality scores)
comp 2 12 0 0 4
  0 cm 22 128
  1 cm 22 128
  2 cm 22 128
  3 mix 14 0 3 12 255
hcomp
  c++ *c=a (store input in M pointed to by C)
  a== 0 if c=0 endif (reset M at newline)
  d=0 b=c hash *d=a a=c a>>= 3 hashd
  d++ a=0 b-- hash *d=a
    b-- a=*b a>>= 5 hashd
  d++ *d=0 b-- a=*b hashd
    b-- a=*b a>>= 4 hashd
  d++ a=*c a>>= 3 *d=0 hashd
    a=c a> 3 if a>>= 5 a+= 4 endif hashd
  halt
post 0 end

Quality scores use a mix of 3 direct context models. The first
uses the previous byte and the column number excluding the
low 3 bits as the context hash. The second model uses the second byte
and the high 3 bits of the third byte back as the context hash.
The third model uses the 4'th byte and the high 4 bits of
the 5'th byte back as context hash. The mixer uses a 14 bit
context consisting of the current partial byte and the column
number with the high 5 bits dropped for column numbers above 3.


ALIGNMENT MODELING

(fxa.cfg to model reference matches)
comp 0 0 0 0 1
  0 cm 20 255
hcomp
  c++ b=a
  a== 0 if a=c a== 1 if c=0 endif endif
  a=c a> 7 if c=0 endif
  a< 6 if
    a=b a>>= 2 a<<= 5 a+=c
  else
    a=c
  endif
  a<<= 9 *d=a
  halt
post 0 end

Reference matches (if present) use a stationary order 0 model with
the parse state (0..7) as context. States 0..3 expect a mismatch
byte and 4..7 expect a pointer byte. States 0..5 also use
the previous byte as context with the low 2 bits discarded.

The ZPAQ archives are each saved as a single segment in a single block
without a locator tag, filename, comment, or checksum. No post-processing
is used. The ZPAQL code used for each of the 4 files is as follows:

Each of the 3 or 4 ZPAQ models is compressed or decompressed in parallel
in separate threads from or to temporary files, which are deleted
when done.

c: input -> output.fx? -> output.fx?.zpaq  (delete output.fx?)
d: input.fx?.zpaq -> input.fx? -> output   (delete input.fx?)
e: input -> output.fx?
f: input.fx? -> output

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <time.h>
#include <pthread.h>
#include "libzpaq.h"
#include "fastq_qscale.h"
using std::string;

const int N=4096; // max FASTQ line length
const int BUCKET=8;

// print error message and exit (may be called by libzpaq)
void libzpaq::error(const char* msg) {
  fprintf(stderr, "fastqz error: %s\n", msg);
  exit(1);
}
using libzpaq::error;

// I/O for libzpaq
struct File: public libzpaq::Reader, public libzpaq::Writer {
  FILE* f;
  int get() {return getc(f);}
  void put(int c) {putc(c, f);}
  int read(char* buf, int n) {return fread(buf, 1, n, f);}
  void write(const char* buf, int n) {fwrite(buf, 1, n, f);}
};

// Thread argument
struct Job {
  int id;  // model 0..2
  string input, output;  // filenames
};

// Thread to compress job.input to job.output using model job.id
void* compress(void *arg) {
  Job& job=*(Job*)arg;
  printf("compressing %s\n", job.input.c_str());

  // Models for fxh, fxb, fxq files
  // Byte codes generated by "zpaq -mfx? l" using fx?.cfg above
  static char hcomp[4][76]={
  {64,0,3,8,0,0,5,2,20,-128,2,22,-128,3,18,3,
  19,7,13,0,4,24,-1,0,104,17,-33,0,47,1,20,28,
  52,74,66,60,68,60,10,68,60,70,25,112,10,68,60,70,
  25,112,10,68,60,70,25,112,10,68,60,25,66,-49,8,112,
  56,0},
  {55,0,3,3,0,0,7,2,9,-1,2,18,-1,2,25,-1,  // fxb
  3,22,8,23,3,4,26,28,7,8,0,6,12,-1,0,17,
  104,74,4,95,1,59,112,10,25,59,112,10,25,59,112,10,
  25,59,112,10,25,59,112,56,0},
  {74,0,2,12,0,0,4,2,22,-128,2,22,-128,2,22,-128,  // fxq
  7,14,0,3,12,-1,0,17,104,-33,0,47,1,20,28,74,
  59,112,66,-41,3,60,25,4,10,59,112,10,68,-41,5,60,
  25,52,10,68,60,10,68,-41,4,60,25,69,-41,3,52,60,
  66,-17,3,47,4,-41,5,-121,4,60,56,0},
  {45,0,0,0,0,0,1,2,20,-1,0,17,72,-33,0,47,
  6,66,-33,1,47,1,20,66,-17,7,47,1,20,-25,6,47,
  8,65,-41,2,-49,5,-126,63,1,66,-49,9,112,56,0}};

  // Compress input to output, then delete input
  libzpaq::Compressor co;
  File in, out;
  in.f=fopen(job.input.c_str(), "rb");
  if (!in.f) perror(job.input.c_str()), exit(1);
  out.f=fopen(job.output.c_str(), "wb");
  if (!out.f) perror(job.output.c_str()), exit(1);
  co.setInput(&in);
  co.setOutput(&out);
  co.startBlock(hcomp[job.id]);
  co.startSegment();
  co.postProcess();
  co.compress();
  co.endSegment();
  co.endBlock();
  fclose(out.f);
  fclose(in.f);
  remove(job.input.c_str());
  printf("compressed %s\n", job.output.c_str());
  return 0;
}

// Thread to decompress job.input to job.output
void* decompress(void *arg) {
  Job& job=*(Job*)arg;
  printf("decompressing %s\n", job.input.c_str());
  File in, out;
  in.f=fopen(job.input.c_str(), "rb");
  if (!in.f) perror(job.input.c_str()), exit(1);
  out.f=fopen(job.output.c_str(), "wb");
  if (!out.f) perror(job.output.c_str()), exit(1);
  libzpaq::decompress(&in, &out);
  fclose(out.f);
  fclose(in.f);
  printf("decompressed %s\n", job.output.c_str());
  return 0;
}

// hash 64 bits to 32 bits
unsigned int hash(unsigned long long hl) {
  return (hl*12345679123456789ull)>>32;
}

// Return the positions of the first 4 mismatches between bbuf[0..len-1]
// and ref[h/4...] (incrementing by dir=(+1,-1)), packed LSB first.
// If there are less than 4 mismatches, use len.
int rmatch(libzpaq::Array<unsigned char>& ref, unsigned int h,
          unsigned char* bbuf, int len, int dir) {
  int i, j, score=0;
  if (len>126) len=126;
  for (i=j=0; i<len && j<4; h+=dir, ++i)
    if (((ref[h/4]>>(6-h%4*2))&3)!=(dir>0?bbuf[i]:3-bbuf[i]))
      score+=i<<(j++*8);
  for (; j<4; ++j)
    score+=len<<(j*8);
  return score;
}

// read reference file into ref
void readref(libzpaq::Array<unsigned char>& ref, const char* filename) {
  FILE* in=fopen(filename, "rb");
  if (!in) perror(filename), exit(1);
  fseek(in, 0, SEEK_END);
  int rlen=ftell(in);
  if (rlen<0 || rlen>=(1<<30))
    error("reference must be smaller than 1 GB");
  rewind(in);
  ref.resize(rlen+N*2);  // pad extra N bytes at each end
  if (int(fread(&ref[N], 1, rlen, in))!=rlen) error("ref read error");
  printf("%s: length=%d bytes\n", filename, rlen);
  fclose(in);
}

void readindex(libzpaq::Array<unsigned int>& index, const char* filename) {
  FILE* in=fopen(filename, "rb");
  if (!in) perror(filename), exit(1);
  int indexSize = (1<<27)+BUCKET;
  index.resize(indexSize);
  if (int(fread(&index[0], sizeof(unsigned int),index.isize(), in))!=indexSize) error("index read error");
  fclose(in);
}

#define openOutputFiles(isref,outPrefix,out) \
do{ \
  int index; \
  for (index=0; index<3+isref; ++index) { \
    string fn=string(outPrefix)+".fx"+"hbqa"[index]; \
    out[index]=fopen(fn.c_str(), "wb"); \
    if (!out[index]) perror(fn.c_str()), exit(1); \
  } \
}while(0)

#define getReadLength(n,in,infile,outfile) \
do{ \
  int i,j,c; \
  for (i=j=n=0; (c=gzgetc(in))!=EOF && !n; ++i) { \
    if (c==10 && j) n=i-j-1; \
    else if (c==10) j=i; \
  } \
  if (n<1 || n>=N) error("read length must be 1..4095"); \
  printf("encoding %s -> %s read length %d\n",infile, outfile, n); \
  gzrewind(in); \
}while(0)

#define saveReadLength(n,out) \
do{ \
  putc(n>>8, out); \
  putc(n&255, out); \
}while(0)

#define encodeHeader(in,hbuf,out) \
do{ \
  int i, j, k, len, c; \
  for (i=j=k=len=0; (c=in[i])!=0; ++i) { \
    if (i>=N) error("Line too long\n"); \
    if (c!=hbuf[i] && isdigit(c) && isdigit(hbuf[i]) && j<254 && i<254 && i==len && (!j || j==i)) { \
      int d=k*10+c-hbuf[i]; \
      if (d>0 && d<254) hbuf[i]=c, k=d, j=i+1; \
    } \
    if (c==hbuf[i] && i==len && len<254) ++len; \
    hbuf[i]=c; \
  } \
  putc(j+(j==0), out); \
  putc(k+1, out); \
  putc(len+1, out); \
  for (j=len; j<i; ++j) putc(hbuf[j], out); \
  putc(0, out); \
}while(0)

#define readBase(in,len,n,hbuf) \
do{ \
  int i, j, c; \
  memset(hbuf,0,N); \
  for (i=0, len=0; (c=in[i])!=0; ++i) { \
    if (c==EOF) error("unexpected EOF"); \
    if (c!='N') { \
      j=(c=='A')+(c=='C')*2+(c=='G')*3+(c=='T')*4; \
      if (!j) error("expected base A,C,G,T,N"); \
      bbuf[len++]=j-1; \
    } \
    hbuf[i]=c; \
  } \
  if (i!=n) error("wrong number of base calls"); \
}while(0)

#define outputMismatch(len,ismatch,bm,bbuf,base,out) \
do{ \
  int i,j; \
  for (i=0; i<len; ++i) { \
    if (!ismatch || i>=(bm>>24&255) || i==(bm>>16&255) || i==(bm>>8&255) || i==(bm&255)) { \
      j="\x01\x03\x04\x02"[bbuf[i]]; \
      if (base*4+j>255) putc(base, out), base=0; \
      base=base*4+j; \
    } \
  } \
}while(0)

#define mapping(len,n,bbuf,index,ref,matches,match_sum,base_sum,base,out1,out2) \
do{ \
  int i, j, k; \
  unsigned long long hl=0; \
  unsigned int bptr=0; \
  int bdir=1,bm=0; \
  bool ismatch=false; \
  for (j=0; j<len; ++j) {  \
    hl=hl*4+bbuf[j]; \
    if (j>=31) { \
      unsigned int h=hash(hl); \
      unsigned int hi=h&0x7ffffff; \
      for (k=0; k<BUCKET && index[hi+k]; ++k) { \
        int m=0; \
        if ((index[hi+k]^h)<0x8000000) { \
          unsigned int ptr=(index[hi+k]&0x7ffffff)*32+31-j; \
          ++matches[n+1]; \
          m=rmatch(ref, ptr, bbuf, len, 1); \
          if (m>bm) bm=m, bptr=ptr; \
        } \
      } \
    } \
  } \
  for (hl=0,j=len-1; j>=0; --j) { \
    hl=hl*4+3-bbuf[j]; \
    if (j<=len-32) { \
      unsigned int h=hash(hl); \
      unsigned int hi=h&0x7ffffff; \
      for (k=0; k<BUCKET && index[hi+k]; ++k) { \
        int m=0; \
        if ((index[hi+k]^h)<0x8000000) { \
          unsigned int ptr=(index[hi+k]&0x7ffffff)*32+31+j; \
          ++matches[n+2]; \
          m=rmatch(ref, ptr, bbuf, len, -1); \
          if (m>bm) bm=m, bptr=ptr, bdir=-1; \
        } \
      } \
    } \
  } \
  ++matches[bm>>24&127]; \
  match_sum+=(bm>>24)&127; \
  match_sum-=(bm^bm<<8)>0xffffff; \
  match_sum-=(bm<<8^bm<<16)>0xffffff; \
  match_sum-=(bm<<16^bm<<24)>0xffffff; \
  base_sum+=len; \
  ismatch=(bm>>23)>=len; \
  if (!ismatch) \
    putc(0, out1); \
  else { \
    putc(1+bm+128*(bdir<0), out1); \
    putc(1+(bm>>8), out1); \
    putc(1+(bm>>16), out1); \
    putc(1+(bm>>24), out1); \
    putc(bptr>>24, out1); \
    putc(bptr>>16, out1); \
    putc(bptr>>8, out1); \
    putc(bptr, out1); \
  } \
  outputMismatch(len,ismatch,bm,bbuf,base,out2); \
}while(0)


#define encodeQual(in,quality,hbuf,n,out) \
do{ \
  int i, j, k, len, c; \
  for (len=0,j=k=0,i=0; (c=in[i])!=0; ++i, k=j, j=c) { \
    if (c!=10 && (c<33 || c>104)) error("expected quality score in 33..104"); \
    if (hbuf[i]=='N') c=33; \
    if (quality>1 && c>35) c-=(c-35)%quality; \
    if (c==35 && (len==0 || j==35)) ++len; \
    else if (len==0 && c>=64 && c<=71) ++len; \
    else if (len==1 && c>=68 && c<=71 && j>=68 && j<=71) ++len; \
    else if (len>=2 && len<55 && k==71 && j==71 && c==71) ++len; \
    else if (c==10 && (len==0 || j==35)) break; \
    else {   \
      ++len;   \
      while (len>1 && j==35) \
        putc(3, out), --len; \
      if (len>3 && j==71 && k==71) \
        putc(199+len, out), len=1; \
      if (len==3) { \
        if (c>=68 && c<=71) \
          putc(137+(k-68)+4*(j-68)+16*(c-68), out), len=0; \
        else \
          putc(73+(k-64)+8*(j-64), out), len=1; \
      } \
      if (len==2) { \
        if (c>=64 && c<=71) putc(73+(j-64)+8*(c-64), out), len=0; \
        else putc(j-32, out), len=1; \
      } \
      if (len==1) { \
        if (c==10) break; \
        if (c!=35 && (c<64 || c>71)) putc(c-32, out), len=0; \
      } \
    } \
  } \
  putc(0, out); \
  if (i!=n) error("wrong number of quality scores"); \
}while(0)

#define print_match_statistics(base_sum,match_sum,matches,n) \
do{ \
  if (base_sum) { \
    printf("matches[0..%d+2]=", n); \
    int i; \
    for (i=0; i<=n+2; ++i) { \
      printf("%d ", matches[i]); \
      if (i%10==0) printf("\n"); \
    } \
    printf("\nMatched %d of %d bases (%1.2f%%)\n",match_sum, base_sum, match_sum*100.0/base_sum); \
  } \
}while(0)

typedef struct _thread_arg_
{
  Fastq *chunk;
  int n;
  int isref;
  int quality;
  const char *outfilePrefix;
  libzpaq::Array<unsigned int> index;
  libzpaq::Array<unsigned char> ref;
  int chunkLen ;
}thread_arg_t;

void *threadEncode64(void *argument){
  thread_arg_t *thread_arg = (thread_arg_t *)argument;
  unsigned char hbuf[N]={0};  // previous header
  unsigned char bbuf[N]={0};  // one sequence
  int bi,len,base=0;
  int matches[N+3]={0},match_sum=0, base_sum=0;
  FILE *out[4];  // fastq, fxh, fxb, fxq, fxa
  openOutputFiles(thread_arg->isref,thread_arg->outfilePrefix,out);
  saveReadLength(thread_arg->n,out[0]);
  char *p=NULL;
  for (bi=0; bi < thread_arg->chunkLen; ++bi){
    for (p=(thread_arg->chunk+bi)->quality;*p!=10;p++) *p-=31;
    encodeHeader((thread_arg->chunk+bi)->name,hbuf,out[0]);
    readBase((thread_arg->chunk+bi)->seq,len,thread_arg->n,hbuf);
    mapping(len,thread_arg->n,bbuf,thread_arg->index,thread_arg->ref,matches,match_sum,base_sum,base,out[3],out[1]);
    encodeQual((thread_arg->chunk+bi)->quality,thread_arg->quality,hbuf,thread_arg->n,out[2]);
    free_Fastq(thread_arg->chunk+bi);
  }
  putc(base, out[1]);
  int i;
  for (i=2+thread_arg->isref; i>=0; --i) fclose(out[i]);
  print_match_statistics(base_sum,match_sum,matches,thread_arg->n);
  pthread_exit(NULL);
}

void *threadEncode33(void *argument){
  thread_arg_t *thread_arg = (thread_arg_t *)argument;
  unsigned char hbuf[N]={0};  // previous header
  unsigned char bbuf[N]={0};  // one sequence
  int bi,len,base=0;
  int matches[N+3]={0},match_sum=0, base_sum=0;
  FILE *out[4];  // fastq, fxh, fxb, fxq, fxa
  openOutputFiles(thread_arg->isref,thread_arg->outfilePrefix,out);
  saveReadLength(thread_arg->n,out[0]);
  for (bi=0; bi < thread_arg->chunkLen; ++bi){
    encodeHeader((thread_arg->chunk+bi)->name,hbuf,out[0]);
    readBase((thread_arg->chunk+bi)->seq,len,thread_arg->n,hbuf);
    mapping(len,thread_arg->n,bbuf,thread_arg->index,thread_arg->ref,matches,match_sum,base_sum,base,out[3],out[1]);
    encodeQual((thread_arg->chunk+bi)->quality,thread_arg->quality,hbuf,thread_arg->n,out[2]);
    free_Fastq(thread_arg->chunk+bi);
  }
  putc(base, out[1]);
  int i;
  for (i=2+thread_arg->isref; i>=0; --i) fclose(out[i]);
  print_match_statistics(base_sum,match_sum,matches,thread_arg->n);
  pthread_exit(NULL);
}

int main(int argc, char** argv) {
  // Start timer
  clock_t start=clock();
  // Check command line: {c|d|e|f} input output
  if (argc<4) {
    printf("fastqz v1.5 FASTQ compressor\n"
    "(C) 2012, Dell Inc. Written by Matt Mahoney. Compiled %s.\n"
    "Licensed under BSD 2 clause license\n"
    "\n"
    "Usage: fastqz command input output [reference]\n"
    "Commands\n"
    "  c[Q] - compress input to output.fx?.zpaq (? = {h,b,q})\n"
    "  d    - decompress input.fx?.zpaq to output\n"
    "  e[Q] - encode (fast) input to output.fx? (? = {h,b,q})\n"
    "  f    - fast decode input.fx? to output\n"
    "Use Q to quantize quality values to steps of size Q for better but\n"
    "lossy compression. Default is c1 or e1 (lossless).\n"
    "Use fapacks to create a reference genome from FASTA files\n",
    __DATE__);
    exit(1);
  }
  const char cmd=argv[1][0]; // c,d,e,f
  int quality=atoi(argv[1]+1);
  if (quality<1) quality=1;
  const int isref=argc>4;    // 1 if a reference file supplied
  int ncpu = getNumCores(), 
  blockSize = argc>6 ? atoi(argv[6]) : 640000,
  threads = argc>5 ? atoi(argv[5]) : ncpu;
//  if (threads > ncpu) threads= ncpu;
  //const int BUCKET=8;        // index bucket size
  libzpaq::Array<unsigned char> ref;  // copy of packed reference genome
  libzpaq::Array<unsigned int> index; // hash table index to ref

  // Encode
  if (cmd=='e' || cmd=='c') {
    // Read reference file
    if (isref) {
      readref(ref, argv[4]);
      printf("load ref at %1.2f seconds\n", double(clock()-start)/CLOCKS_PER_SEC);
      string fn=string(argv[4])+".index";
      readindex(index, fn.c_str());
      printf("load index at %1.2f seconds\n", double(clock()-start)/CLOCKS_PER_SEC);
    }

    int n, i, j, k;

    gzFile in=open_input_stream(argv[2]);
    if (!in) perror(argv[2]), exit(1);
    int qscale_index = fastq_qualscale(argv[2], in);
    string fn=string(argv[3])+".fxs";
    FILE *out_qscale=fopen(fn.c_str(), "wb");
    if (!out_qscale) perror(fn.c_str()), exit(1);
    fprintf(out_qscale,"%d",qscale_index);
    fclose(out_qscale);
    printf("check quality scale at %1.2f seconds\n", double(clock()-start)/CLOCKS_PER_SEC);

    // Save read length, n
    getReadLength(n,in,argv[2],argv[3]);

    int line=0,eofFlag=0,blockCount=0;
    char *LineBuf = (char *)calloc(N,sizeof(char)),ret[32];
    Fastq **chunk=(Fastq **)calloc(threads,sizeof(Fastq *));
    char **retBuf = (char **)calloc(threads,sizeof(char *));
    int *chunkLen = (int *)calloc(threads,sizeof(int));
    for (i=0;i<threads;i++){
      retBuf[i]=(char *)calloc(512,sizeof(char));
      chunk[i]=(Fastq *)calloc(blockSize,sizeof(Fastq));
    }

    for (line=0,k=0; 1; line++) {
      int blockIndex = line % blockSize;
      if ((!blockIndex && line) || eofFlag) {
        chunkLen[k++ % threads]=!eofFlag ? blockSize : blockIndex-1;
        if (!eofFlag||(eofFlag && blockIndex>1)) ++blockCount;
        int chunkIndex = blockCount % threads;
        if ((!chunkIndex && line) || eofFlag) {
          int chunkCount = !eofFlag ? threads : chunkIndex;
          pthread_t *thread_handle = (pthread_t *)malloc(chunkCount * sizeof(pthread_t));
          thread_arg_t *thread_args = (thread_arg_t *)malloc(chunkCount* sizeof(thread_arg_t));
          for (j=0;j<chunkCount ;j++) {
            sprintf(retBuf[j],"%s.%d",argv[3],((blockCount-1)/threads)*threads+j);
            (thread_args+j)->chunk=chunk[j];
            (thread_args+j)->n=n;
            (thread_args+j)->isref=isref;
            (thread_args+j)->quality=quality;
            (thread_args+j)->outfilePrefix = retBuf[j];
            (thread_args+j)->index = index;
            (thread_args+j)->ref = ref;
            (thread_args+j)->chunkLen = chunkLen[j];
            pthread_create(thread_handle+j,NULL,(qscale_index >1 || qscale_index<0) ? threadEncode64 : threadEncode33,thread_args+j);
          }
          for(j = 0;j < chunkCount;j++){
            pthread_join(*(thread_handle + j),NULL);
          }
          free(thread_args);
          free(thread_handle);
        }
      }
      if (eofFlag) break;
      eofFlag = readNextNode(in,LineBuf,&chunk[blockCount % threads][blockIndex]);
    }
    gzclose(in);
    index.resize(0);
    ref.resize(0);
    printf("encode file at %1.2f seconds\n", double(clock()-start)/CLOCKS_PER_SEC);

    // compress each temporary file to .zpaq in a separate thread
    
    if (cmd=='c') {
      for (k=0;k<blockCount;k++){
        pthread_t tid[4];
        pthread_attr_t attr; // thread joinable attribute
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        Job job[4];
        for (i=0; i<3+isref; ++i) {
          job[i].id=i;
          int2str(k,10,ret);
          job[i].input=string(argv[3])+'.'+string(ret)+".fx"+"hbqa"[i];
          job[i].output=job[i].input+".zpaq";
          pthread_create(&tid[i], &attr, compress, (void*)&job[i]);
        }

        // wait until all jobs are done
        for (i=0; i<3+isref; ++i) {
          void* status;
          pthread_join(tid[i], &status);
        }
      }
    }
  }

  // decode
  else if (cmd=='d' || cmd=='f') {

    // decompress .zpaq
    Job job[4];
    if (cmd=='d') {
      pthread_t tid[4];
      pthread_attr_t attr; // thread joinable attribute
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      for (int i=0; i<3+isref; ++i) {
        job[i].id=i;
        job[i].output=string(argv[2])+".fx"+"hbqa"[i];
        job[i].input=job[i].output+".zpaq";
        pthread_create(&tid[i], &attr, decompress, (void*)&job[i]);
      }

      // wait until all threads are done
      for (int i=0; i<3+isref; ++i) {
        void* status;
        pthread_join(tid[i], &status);
      }
    }

    // read reference
    if (isref) readref(ref, argv[4]);

    // open  files
    FILE *in[5], *out;  // fxh, fxb, fxq, fxa, fastq
    int i, j, k, c, n;
    for (i=0; i<3+isref; ++i) {
      string fn=string(argv[2])+".fx"+"hbqa"[i];
      in[i]=fopen(fn.c_str(), "rb");
      if (!in[i]) perror(fn.c_str()), exit(1);
    }
    out=fopen(argv[3], "wb");
    if (!out) perror(argv[3]), exit(1);

    string fn=string(argv[2])+".fxs";
    in[4]=fopen(fn.c_str(), "rb");
    if (!in[4]) perror(fn.c_str()), exit(1);
    int qscale_index=0;
    int matched = fscanf(in[4],"%d",&qscale_index);
    if (matched ) fprintf(stderr,"qscale_index: %d\n",qscale_index);
    fclose(in[4]);
    remove(fn.c_str());

    // get read length, n
    n=getc(in[0]);
    n=n*256+getc(in[0]);
    printf("decoding %s -> %s read length %d\n",
      argv[2], argv[3], n);
    if (n<1 || n>=N) error("bad read length");

    // decode
    int base=0;
    unsigned char hbuf[N]={0}, qbuf[N]={0};
    while (1) {

      // decode header
      j=getc(in[0])-1;  // index of last digit of number to adjust
      if (j==EOF-1) break;
      k=getc(in[0])-1;  // amount to add
      i=getc(in[0])-1;  // number of matched bytes after adjustment
      if (j<0 || k<0 || i<0) error("bad header");
      for (; i<N && (c=getc(in[0]))!=EOF && c; ++i) hbuf[i]=c;
      for (; k && j>=0; --j, k/=10) {
        int d=k%10;
        hbuf[j]+=d, k-=d;
        if (hbuf[j]>'9') hbuf[j]-=10, k+=10;
      }
      for (j=0; j<i; ++j) putc(hbuf[j], out);
      putc(10, out);

      // read quality scores and save in qbuf[0..n-1]
      // 0 -> pad with 35 and end
      // c=1..72 -> c+32
      // c=73..136 -> (c-73)%8+64, (c-73)/8+64
      // c=137..200 -> (c-137)%4+68, (c-137)/4%4+68, (c-137)%16+68
      // c=201..255 -> 71 repeated c-200 times
      for (i=0;;) {
        c=getc(in[2]);
        if (c==EOF) error("unexpected end of .fxq");
        if (i>n) error("missing .fxq terminator");
        if (c==0) { // end of line
          for (; i<n; ++i) qbuf[i]=35;
          break;
        }
        else if (c>=201 && i+c-200<=n)
          while (c-->200) qbuf[i++]=71;
        else if (c>=137 && c<=200 && i<n-2) {
          c-=137;
          qbuf[i++]=(c&3)+68;
          qbuf[i++]=((c>>2)&3)+68;
          qbuf[i++]=((c>>4)&3)+68;
        }
        else if (c>=73 && c<=136 && i<n-1) {
          c-=73;
          qbuf[i++]=(c&7)+64;
          qbuf[i++]=((c>>3)&7)+64;
        }
        else if (c>=1 && c<=72 && i<n) {
          qbuf[i++]=c+32;
        }
        else error (".fxq code overflow");
      }
      if (i!=n) error("incorrect .fxq read length");

      // decode match to reference
      unsigned int bptr=0;  // pointer to match in ref
      int bdir=0;  // read direction
      int miss1=0, miss2=0, miss3=0, miss4=0;  // mismatches, ascending order
      if (isref) {
        miss1=getc(in[3]);
        if (miss1==EOF) error("unexpcted EOF in .fxa");
        if (miss1) {
          if (miss1>=128) miss1-=128, bdir=-1;
          else bdir=1;
          --miss1;
          miss2=getc(in[3])-1;
          miss3=getc(in[3])-1;
          miss4=getc(in[3])-1;
          bptr=getc(in[3]);
          bptr=bptr*256+getc(in[3]);
          bptr=bptr*256+getc(in[3]);
          bptr=bptr*256+getc(in[3]);
        }
      }

      // decode bases
      for (i=k=0; i<n; ++i) {
        if (qbuf[i]==33)
          putc('N', out);
        else if (bdir && k!=miss1 && k!=miss2 && k!=miss3 && k<miss4) {
          if (bptr/4>=ref.size()) error(".fxa pointer out of bounds");
          j=(ref[bptr/4]>>(6-bptr%4*2))&3;
          bptr+=bdir;
          if (bdir<0) j=3-j;
          putc("ACGT"[j], out);
          ++k;
        }
        else {
          while (base==0) {
            base=getc(in[1]);
            if (base==EOF) error("unexpected end of .fxb");
          }
          if (base>84) j=(base-21)>>6, base-=j*64;
          else if (base>20) j=(base-5)>>4, base-=j*16;
          else if (base>4) j=(base-1)>>2, base-=j*4;
          else j=base, base=0;
          putc(" ATCG"[j], out);
          ++k;
          bptr+=bdir;
        }
      }
      putc(10, out);

      // write empty second header
      putc('+', out);
      putc(10, out);

      // write quality scores
      if (qscale_index >1 || qscale_index <0){
        for (i=0; i<n; ++i){
          if (qbuf[i]==33) qbuf[i]=35;
          putc(qbuf[i]+31, out);
        }
      }else{
        for (i=0; i<n; ++i){
          if (qbuf[i]==33) qbuf[i]=35;
          putc(qbuf[i], out);
        }
      }
      putc(10, out);
    }
    fclose(out);
    for (i=2+isref; i>=0; --i) fclose(in[i]);

    // delete temporary files
    if (cmd=='d')
      for (int i=0; i<3+isref; ++i)
        remove(job[i].output.c_str());

    // show results
    printf("decoded %s\n", argv[3]);
  }
  printf("%1.2f seconds\n", double(clock()-start)/CLOCKS_PER_SEC);
  return 0;
}
