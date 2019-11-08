/*************************************************************************
 * C <cstdlib> library qsort function.
*************************************************************************/
template <typename T>
T cmp(const void *a, const void *b) { return (*(T*)a - *(T*)b); }

template <typename T>
void qSort(T* begin, T* end) { qsort(begin, std::distance(begin, end), sizeof(T), cmp); }

/*************************************************************************
 * C++ STL <algorithm> library qsort function.
*************************************************************************/
template<typename T>
T compare(const void *a, const void *b) { return (*(T*)a - *(T*)b); }

template<typename T>
void stlQSort(T* begin, T* end) { std::qsort(begin, std::distance(begin, end), sizeof(T), compare); }
