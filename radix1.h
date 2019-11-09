const int bits = 31;
const int radix = 10;
const int digits = 10;

template<class T>
void radixSort(T data[], const int n) 
{
  register int d, j, k, factor;
  Queue<T> queues[radix];

  for (d = 0, factor = 1; d < digits; factor *= radix, d++) 
  {
    for (j = 0; j < n; j++)
      queues[(data[j] / factor) % radix].enqueue(data[j]);

    for (j = k = 0; j < radix; j++)
      while (!queues[j].empty())
        data[k++] = queues[j].dequeue();
  }
}

void bitRadixSort(long data[], const int n, int b) 
{
  int pow2b = 1;
  pow2b <<= b;
  int i, j, k, pos = 0, mask = pow2b - 1;
  int last = (bits % b == 0) ? (bits / b) : (bits / b + 1);
  Queue<long> *queues = new Queue<long>[pow2b];

  for (i = 0; i < last; i++) 
  {
    for (j = 0; j < n; j++)
      queues[(data[j] & mask) >> pos].enqueue(data[j]);
    mask <<= b;
    pos = pos + b;

    for (j = k = 0; j < pow2b; j++)
      while (!queues[j].empty())
        data[k++] = queues[j].dequeue();
  }
}
