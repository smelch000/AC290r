#include <iostream>
#include <algorithm>
#include "double3.h"

std::istream & operator>> (std::istream & in, double3 & val) {
  in >> val.x >> val.y >> val.z;
  return in;
}

std::ostream & operator<< (std::ostream & out, double3 const & val) {
  out << val.x << " " << val.y << " " <<  val.z;
  return out;
}

void double3::from_string(std::string const & s, int skipwords, bool replaceDecimal) {
  std::string ss;
  ss=s;
  if (replaceDecimal) {
    std::replace(ss.begin(),ss.end(),',','.');
  }
  std::istringstream iss(ss);
  std::string sk;
  for (int i=0;i<skipwords; i++) iss >> sk;
  iss >> *this;
}

std::string double3::to_string() {
  std::ostringstream iss;
  iss << *this;
  return (iss.str());
}
