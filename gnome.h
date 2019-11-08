/*
// gnome sort pseudocode
function gnomeSort(a[0..size-1])
  i := 1
  j := 2
  while i < size do
    if a[i-1] <= a[i] then
      // for descending sort, use >= for comparison
      i := j
      j := j + 1
    else
      swap a[i-1] and a[i]
      i := i - 1
      if i = 0 then
        i := j
        j := j + 1
      endif
    endif
  done

template<typename RandomAccessIterator>
void gnome_sort(RandomAccessIterator begin, RandomAccessIterator end) 
{
  auto i = begin + 1;
  auto j = begin + 2;

  while (i < end) {
  if (!(*i < *(i - 1))) {
    i = j;
    ++j;
  } else {
    std::iter_swap(i - 1, i);
    --i;
    if (i == begin) {
    i = j;
    ++j;
    }
  }
  }
}
*/

template <typename T>
void gnomeSort(T data[], std::size_t size)
{
  T temp;

  for (std::size_t i = 0; i < size - 1; i++)
  {
    temp = data[i];

    if (data[i] > data[i + 1])
    {
      data[i] = data[i + 1];
      data[i + 1] = temp;

      for (std::size_t j = i; j > 0; j--)
      {
        temp = data[j];

        if (data[j] < data[j - 1])
        {
          data[j] = data[j - 1];
          data[j - 1] = temp;
        }
      }
    }
  }
}
