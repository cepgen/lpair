#ifndef _UTILS_H
#define _UTILS_H

#include <ctime>

/**
 * An object which enables to extract the processing time between two steps in
 * this software's flow
 */
class Timer
{
 public:
  inline Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }
  /**
   * Get the time elapsed since the last @a reset call (or class construction)
   * @return The elapsed time in seconds
   */
  inline double elapsed() {
    clock_gettime(CLOCK_REALTIME, &end_);
    return end_.tv_sec -beg_.tv_sec+(end_.tv_nsec - beg_.tv_nsec)/1000000000.;
  }
  /**
   * @brief Resets the clock counter
   */
  inline void reset() {
    clock_gettime(CLOCK_REALTIME, &beg_);
  }
 private:
  /** @brief Timestamp marking the beginning of the counter */
  timespec beg_;
  /** @brief Timestamp marking the end of the counter */
  timespec end_;
};

#endif
