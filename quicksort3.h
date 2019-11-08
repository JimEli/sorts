/*************************************************************************
 * Implements a basic quick sort algorithm.
 *************************************************************************/
template<class T>
void quick_sort(T* begin, int first, int last)
{
  int lower = first + 1, upper = last;

  //std::swap(*(begin + first), *(begin + (first + last)/2));
  T bound = *(begin + first);

  while (lower <= upper)
  {
    while (*(begin + lower) < bound)
      lower++;

    while (bound < *(begin + upper))
      upper--;

    //if (lower < upper)
      //std::swap(*(begin + lower++), *(begin + upper--));
    //else
      //lower++;
    if (lower < upper)
    {
      std::swap(*(begin + lower), *(begin + upper));
      upper--;
    }
    lower++;
  }

  std::swap(*(begin + upper), *(begin + first));

  if (first < upper - 1)
    quick_sort(begin, first, upper - 1);

  if (upper + 1 < last)
    quick_sort(begin, upper + 1, last);
}

template<class T>
void quickSort(T* begin, T* end)
{
  const int n = std::distance(begin, end);

  if (n < 2)
    return;

  int max = 0;
  // Find largest element and put it at the end.
  for (int i = 1; i < n; i++)
    if (*(begin + max) < *(begin + i))
      max = i;

  // Largest element is now in its final position.
  std::swap(*(begin + n - 1), *(begin + max));
  quick_sort(begin, 0, n - 2);
}
