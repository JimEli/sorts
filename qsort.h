/*************************************************************************
 * C <cstdlib> library qsort function.
*************************************************************************/
template <typename T>
T cmp(const void *a, const void *b) { return (*(T*)a - *(T*)b); }

template <typename T>
void qSort(T* begin, T* end) { qsort(begin, std::distance(begin, end), sizeof(T), cmp); }
