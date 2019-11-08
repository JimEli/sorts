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
