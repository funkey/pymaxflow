#include <Graph.h>

typedef Graph<int>   GraphInt;
typedef Graph<float> GraphFloat;

template <> Graph<int>;
template <> Graph<float>;
