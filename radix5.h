/*************************************************************************
 * Implements least significant digit radix sort algorithm.
 *************************************************************************/
 // Radix sort comparator for 32-bit two's complement integers.
class radixTest 
{
  // Bit position to test.
  const int bit;

public:
  radixTest(int offset) : bit(offset) { }

  bool operator() (int value) const
  {
    // Sign bit?
    if (bit == 31)
      // Negative value sorted to left partition.
      return value < 0;
    else
      // 0 bit sorted to left partition.
      return !(value & (1 << bit));
  }
};

// Least significant digit radix sort.
template <typename T>
void radixSort(T* begin, T* end)
{
  for (int lsb = 0; lsb < 32; ++lsb)
    std::stable_partition(begin, end, radixTest(lsb));
}
