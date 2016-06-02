/**
 * Utility functions for measuring memory usage
 * See discussions at: 
 * http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
 * http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system
 */

#ifndef MEMUTILS_H
#define	MEMUTILS_H


#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

#ifdef __linux__
/**
 * Uses a linux-specific mechanism to determine current process memory usage
 * @return Current used memory (resident set size) in KB
 */
double process_mem_usage()
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   return rss * page_size_kb;
}
#else
double process_mem_usage()
{
    return 0.0;
}
#endif


#endif	/* MEMUTILS_H */

