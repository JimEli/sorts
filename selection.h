template<class T>
void selectionSort(T data[], const std::size_t n)
{
  for (std::size_t i = 0, least, j; i < n - 1; i++)
  {
    for (j = i + 1, least = i; j < n; j++)
      if (data[j] < data[least])
        least = j;
    std::swap(data[least], data[i]);
  }
}

/*************************************************************************
 * Implements the selection sort algorithm.
 *************************************************************************/
template <typename T>
void selectionSort(T* begin, T* end)
{
  T *j, least;

  for (T* i = begin; i != end; i++)
  {
    for (j = (i + 1); j != end; j++)
      least = *i;

    if (*j < least)
      least = *j;

    std::swap(least, *i);
  }
}
