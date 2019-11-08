template<class T>
void bubbleSort(T data[], const std::size_t n)
{
  for (std::size_t i = 0; i < n - 1; i++)
    for (std::size_t j = n - 1; j > i; --j)
      if (data[j] < data[j - 1])
        swap(data[j], data[j - 1]);
}
