template<class T>
void quickSort2(T data[], int first, int last) 
{
  if (last - first < 30)
    insertionSort(data, first, last);
  else 
  {
    int lower = first + 1, upper = last;
    swap(data[first], data[(first + last) / 2]);
    T bound = data[first];

    while (lower <= upper) {
      while (data[lower] < bound)
        lower++;

      while (bound < data[upper])
        upper--;

      if (lower < upper)
        swap(data[lower++], data[upper--]);
      else
        lower++;
    }
    swap(data[upper], data[first]);

    if (first < upper - 1)
      quickSort2(data, first, upper - 1);

    if (upper + 1 < last)
      quickSort2(data, upper + 1, last);
  }
}

template<class T>
void quickSort2(T data[], const int n) 
{
  int i, max;

  if (n < 2)
    return;

  for (i = 1, max = 0; i < n; i++)// find the largest
    if (data[max] < data[i])    // element and put it
      max = i;                // at the end of data[];

  swap(data[n - 1], data[max]);   // largest el is now in its
  quickSort2(data, 0, n - 2);     // final position;
}
