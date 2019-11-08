template<class T>
void insertionSort(T data[], const std::size_t n)
{
  for (std::size_t i = 1, j; i < n; i++)
  {
    T tmp = data[i];

    for (j = i; j > 0 && tmp < data[j - 1]; j--)
      data[j] = data[j - 1];
    data[j] = tmp;
  }
}

/*************************************************************************
 * Implements the insertion sort algorithm.
 *************************************************************************/
template<class T>
void insertionSort(T* begin, T* end)
{
  for (T* i = begin; i != end; i++)
  {
    T tmp = *i, *j;

    for (j = i; j >= begin && tmp < *(j - 1); j--)
      *j = *(j - 1);
    *j = tmp;
  }
}

/*
// Binary insertion sort O(n^2) time O(1) space, binary search to insert instead of swapping
template <typename T>
void insertion_sort(T* begin, T* end, T* cur) 
{
  for (; cur < end; ++cur)
  {
    T pivot = std::move(*cur);

    // upper_bound is binary search for correct insertion position
    const T pos = std::upper_bound(begin, cur, pivot);
    for (T* ins = cur; ins != pos; --ins)
      *ins = std::move(*(ins - 1));
    *pos = std::move(pivot);
  }
}

template <typename T>
void insertionSort(T* begin, T* end) { insertion_sort(begin, end, begin); }
*/
