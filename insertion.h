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
