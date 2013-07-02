/*
 * Debug.h
 *
 *  Created on: 2013-07-02
 *      Author: Michael
 */

#ifndef DEBUG_H_
#define DEBUG_H_
#include <ostream>

class Debug {
 public:

  // Standard stuff
  template<typename T>
  Debug& operator <<(const T &val) {
#ifdef DEBUG
    std::cerr << val;
#endif
    return *this;
  }

  // For things like endl
  Debug& operator <<(std::ostream& (*pf)(std::ostream&)) {
#ifdef DEBUG
    std::cerr << pf;
#endif
    return *this;
  }

  static Debug& getInstance() {
    return debugInstance_;
  }

 private:
  static Debug debugInstance_;
  Debug();
};

#endif /* DEBUG_H_ */
