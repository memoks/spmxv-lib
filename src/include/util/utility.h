
#ifndef UTILITY_H_
#define UTILITY_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"

extern void accumulate(int* arr_inout, int length);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* UTILITY_H_ */
