
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

void accumulate(int* arr_inout, int length)
{
	if(arr_inout == NULL)
		return;

	int i;
	for(i = 1; i < length; ++i)
		arr_inout[i] += arr_inout[i - 1];
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
