template<class T>
void merge(T array1[], T temp[], int first, int last) 
{
  int mid = (first + last) / 2;
  int i1 = 0, i2 = first, i3 = mid + 1;

  while (i2 <= mid && i3 <= last)
    if (array1[i2] < array1[i3])
      temp[i1++] = array1[i2++];
    else
      temp[i1++] = array1[i3++];

  while (i2 <= mid)
    temp[i1++] = array1[i2++];

  while (i3 <= last)
    temp[i1++] = array1[i3++];

  for (i1 = 0, i2 = first; i2 <= last; array1[i2++] = temp[i1++]);
}

template<class T>
void mergeSort(T data[], T temp[], int first, int last) 
{
  if (first < last) 
  {
    int mid = (first + last) / 2;

    mergeSort(data, temp, first, mid);
    mergeSort(data, temp, mid + 1, last);
    merge(data, temp, first, last);
  }
}

template<class T>
void mergeSort(T data[], const int n) 
{
  T *temp = new T[n];
  mergeSort(data, temp, 0, n - 1);
}

/*************************************************************************
 * Implements merge sort algorithm using STL comparison function std::less.
 *************************************************************************/
template <typename T>
void mergeSort(T* begin, T* end)
{
  if (end - begin > 1)
  {
    T* mid = begin + (end - begin) / 2;
    mergeSort(begin, mid);
    mergeSort(mid, end);
    std::inplace_merge(begin, mid, end, std::less<T>());
  }
}
