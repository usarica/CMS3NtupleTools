#ifndef PARTICLEOBJECTHELPERS_H
#define PARTICLEOBJECTHELPERS_H

#include <vector>
#include <algorithm>
#include "ParticleObject.h"


namespace ParticleObjectHelpers{
  template<typename T> bool objHasGreaterPt(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterPt(T const* earlier, T const* later);

  template<typename T> bool objHasGreaterScalarSumPt(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterScalarSumPt(T const* earlier, T const* later);

  template<typename T> void sortByGreaterPt(std::vector<T>& vec);
  template<typename T> void sortByGreaterPt(std::vector<T*>& vec);

  template<typename T> void sortByGreaterScalarSumPt(std::vector<T>& vec);
  template<typename T> void sortByGreaterScalarSumPt(std::vector<T*>& vec);

}

template<typename T> bool ParticleObjectHelpers::objHasGreaterPt(T const& earlier, T const& later){ return (earlier.pt() > later.pt()); }
template<typename T> bool ParticleObjectHelpers::ptrHasGreaterPt(T const* earlier, T const* later){ return (earlier && ((later && earlier->pt() > later->pt()) || !later)); }

template<typename T> bool ParticleObjectHelpers::objHasGreaterScalarSumPt(T const& earlier, T const& later){
  std::vector<ParticleObject const*> deepDaus_earlier; earlier.getDeepDaughters(deepDaus_earlier);
  std::vector<ParticleObject const*> deepDaus_later; later.getDeepDaughters(deepDaus_later);
  float scsumpt_earlier=0; for (auto dau:deepDaus_earlier) scsumpt_earlier += dau->pt();
  float scsumpt_later=0; for (auto dau:deepDaus_later) scsumpt_later += dau->pt();
  return (scsumpt_earlier > scsumpt_later);
}
template<typename T> bool ParticleObjectHelpers::ptrHasGreaterScalarSumPt(T const* earlier, T const* later){
  return (earlier && ((later && ParticleObjectHelpers::objHasGreaterScalarSumPt(*earlier, *later)) || !later));
}

template<typename T> void ParticleObjectHelpers::sortByGreaterPt(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::objHasGreaterPt<T>); }
template<typename T> void ParticleObjectHelpers::sortByGreaterPt(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::ptrHasGreaterPt<T>); }

template<typename T> void ParticleObjectHelpers::sortByGreaterScalarSumPt(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::objHasGreaterScalarSumPt<T>); }
template<typename T> void ParticleObjectHelpers::sortByGreaterScalarSumPt(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::ptrHasGreaterScalarSumPt<T>); }


#endif
