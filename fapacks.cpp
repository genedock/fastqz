/* fapacks.cpp - pack FASTA 4 bases per byte
   includes lowercase a,c,g,t

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

This program produces packed DNA sequences from FASTA files.
The output may be used as a reference genome for the program fastqz.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <time.h>
#include "libzpaq.h"
using std::string;

const int N=4096; // max FASTQ line length
const int BUCKET=8;

// print error message and exit (may be called by libzpaq)
void libzpaq::error(const char* msg) {
  fprintf(stderr, "fastqz error: %s\n", msg);
  exit(1);
}

using libzpaq::error;

// hash 64 bits to 32 bits
unsigned int hash(unsigned long long hl) {
  return (hl*12345679123456789ull)>>32;
}

void writeindex(libzpaq::Array<unsigned int>& index, const char* filename) {
  FILE* out=fopen(filename, "wb");
  if (!out) perror(filename), exit(1);
  fwrite(&index[0], sizeof(unsigned int),index.isize(), out);
  fclose(out);
}

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
  fprintf(stderr,"%s: length=%d bytes\n", filename, rlen);
  fclose(in);
}

void buildindex(libzpaq::Array<unsigned char>& ref,libzpaq::Array<unsigned int>& index, const char* filename) {
  readref(ref, filename);  // read into ref

  // Create an index. Divide ref into groups of 32 bases (8 bytes)
  // and compute a 32 bit hash, h. Use the low 27 bits as a hash index
  // and high 5 bits as a hash checksum. Store the checksum and a
  // 27 bit pointer into ref packed into index[h].
  index.resize((1<<27)+BUCKET);
  int collisions=0;
  for (int i=N; i<=int(ref.size())-N-8; i+=8) {
    unsigned long long hl=0;
    for (int j=0; j<8; ++j) hl=hl<<8|ref[i+j];
    unsigned int h=hash(hl);
    unsigned int hi=h&0x7ffffff;
    int j;
    for (j=0; j<BUCKET && index[hi+j]; ++j);
    if (j==BUCKET) ++collisions;
    else index[hi+j]=(h&0xf8000000)+(i>>3);
  }
  fprintf(stderr,"indexed %s: %d of %lu collisions\n",filename, collisions, ref.size()/8);
}


int main(int argc, char** argv) {
  if (argc<3)
    fprintf(stderr,"To pack FASTA files: fapack output *.fa\n"), exit(1);
  FILE *out=fopen(argv[1], "wb");
  int b=1, c;
  for (int i=2; i<argc; ++i) {
    fprintf(stderr,"%s\n", argv[i]);
    FILE *in=fopen(argv[i], "rb");
    if (!in) continue;
    bool dna=true;
    while ((c=getc(in))!=EOF) {
      if (c=='>') dna=false;
      else if (c==10) dna=true;
      if (islower(c)) c=toupper(c);
      if (dna) {
        if (c=='A') b=b*4;
        if (c=='C') b=b*4+1;
        if (c=='G') b=b*4+2;
        if (c=='T') b=b*4+3;
        if (b>=256) putc(b&255, out), b=1;
      }
    }
    if (in) fclose(in);
  }
  fclose(out);

  libzpaq::Array<unsigned char> ref;  // copy of packed reference genome
  libzpaq::Array<unsigned int> index; // hash table index to ref
  
  buildindex(ref,index,argv[1]);
  string fn=string(argv[1])+".index";
  writeindex(index,fn.c_str());
  return 0;
}


