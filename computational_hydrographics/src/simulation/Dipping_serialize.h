#ifndef CH_DIPPING_SERIALIZE
#define CH_DIPPING_SERIALIZE

#include "Dipping.h"
#include <igl/serialize.h>

void dipping_serialize(const std::string filename, const Dipping& c);
void dipping_deserialize(const std::string filename, Dipping& c);

#endif
