/*************************************************************************
 * C++ STL <algorithm> library sort function. sort() uses an introsort
 * (quicksort which switches to heapsort when the recursion reaches a
 * certain depth).
*************************************************************************/
template<typename T>
void stlSort(T* begin, T* end) { std::sort(begin, end); }
