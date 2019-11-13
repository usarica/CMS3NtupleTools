#ifndef METFILTERINFO_H
#define METFILTERINFO_H

#include <string>
#include <unordered_map>


struct METFilterInfo{
  std::unordered_map<std::string, bool> flag_accept_map;
};


#endif
