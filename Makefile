default: fastqz fapack fapacks

fastqz: fastqz15.cpp libzpaq.3.pod libzpaq.cpp libzpaq.h
	g++ -O2 -msse2 fastqz15.cpp libzpaq.cpp -o $@ -lz -lpthread

fastqz_thread: fastqz15_thread.cpp libzpaq.3.pod libzpaq.cpp libzpaq.h
	g++ -O2 -msse2 fastqz15_thread.cpp libzpaq.cpp -o $@ -lz -lpthread

fastqz_thread_Cache: fastqz15_thread_Cache.cpp libzpaq.3.pod libzpaq.cpp libzpaq.h
	g++ -O2 -msse2 fastqz15_thread_Cache.cpp libzpaq.cpp -o $@ -lz -lpthread

fastqz_thread_CacheFile: fastqz15_thread_CacheFile.cpp libzpaq.3.pod libzpaq.cpp libzpaq.h
	g++ -O2 -msse2 fastqz15_thread_CacheFile.cpp libzpaq.cpp -o $@ -lz -lpthread

fapacks: fapacks.cpp libzpaq.cpp libzpaq.h
	g++ -O2 -msse2 fapacks.cpp libzpaq.cpp -o $@

clean:
	- rm -f fastqz fapack fapacks fastqz_thread fastqz_thread_Cache
