void countingSort(long data[], const long n) 
{
  long i;
  long largest = data[0];
  long *tmp = new long[n];

  for (i = 1; i < n; i++)        // find the largest number
    if (largest < data[i])     // in data and create the array
      largest = data[i];     // of counters accordingly;

  unsigned long *count = new unsigned long[largest + 1];

  for (i = 0; i <= largest; i++)
    count[i] = 0;

  for (i = 0; i < n; i++)        // count numbers in data[];
    count[data[i]]++;

  for (i = 1; i <= largest; i++) // count numbers <= i; 
    count[i] = count[i - 1] + count[i];

  for (i = n - 1; i >= 0; i--)   // put numbers in order in tmp[];
  {
    tmp[count[data[i]] - 1] = data[i];
    count[data[i]]--;
  }

  for (i = 0; i < n; i++)        // transfer numbers from tmp[]
    data[i] = tmp[i];          // to the original array;
}

/*************************************************************************
 * Implements a counting sort algorithm.
 *************************************************************************/
template<typename T>
void countingSort(T* begin, T* end)
{
  T min = *std::min_element(begin, end);
  T max = *std::max_element(begin, end);
  std::vector<unsigned> count((max - min) + 1, 0);

  for (T* i = begin; i != end; ++i)
    ++count[*i - min];

  T* start = begin;

  for (T i = min; i <= max; ++i)
    for (unsigned j = 0; j < count[i - min]; ++j)
      *start++ = i;
}
