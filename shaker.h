/*************************************************************************
 * Implements the cocktail shaker sort algorithm.
 *************************************************************************/
template <typename T>
void shakerSort(T* begin, T* end)
{
  bool swapped = true;
  T* last = end;

  while (begin != last-- && swapped)
  {
    swapped = false;

    for (T* i = begin; i != last; ++i)
    {
      if (*(i + 1) < *i)
      {
        std::swap(*i, *(i + 1));
        swapped = true;
      }
    }

    if (!swapped)
      break;

    swapped = false;

    for (T* i = last - 1; i != begin; --i)
    {
      if (*i < *(i - 1))
      {
        std::swap(*i, *(i - 1));
        swapped = true;
      }
    }

    ++begin;
  }
}
