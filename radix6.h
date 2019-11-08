/*************************************************************************
 * Implements a radix sort algorithm.
 *************************************************************************/
template <typename T>
void radixSort(T* begin, T* end)
{
  T* temp, *first = begin;
  T* buffer = new T[std::distance(begin, end)], *pBuffer = buffer;
  constexpr int DIGIT_BITS = 8;                   // Bits per digit.
  constexpr size_t NUM_BUCKETS = 1 << DIGIT_BITS; // Number of buckets.
  constexpr T BUCKET_MASK = NUM_BUCKETS - 1;      // Mask all bits.
  uint32_t bucketSize[NUM_BUCKETS];
  uint32_t bucketPos[NUM_BUCKETS];

  for (uint32_t j = 0; j < (8 * sizeof(T)); j += DIGIT_BITS)
  {
    std::memset(bucketSize, 0, sizeof(bucketSize));

    for (temp = first; temp != end; ++temp)
      ++bucketSize[(*temp >> j) & BUCKET_MASK];

    if (bucketSize[((*first >> j) & BUCKET_MASK)] != uint32_t(end - first))
    {
      if (j == (8 * sizeof(T) - DIGIT_BITS))
      {
        bucketPos[0] = 0;

        for (uint32_t i = 0; i < NUM_BUCKETS - 1; i++)
          bucketPos[i + 1] = bucketPos[i] + bucketSize[i];

        for (temp = first; temp != end; ++temp)
        {
          const size_t digit = (*temp >> j) & BUCKET_MASK;
          *(pBuffer + bucketPos[digit]++) = *temp;
        }
      }
      else
      {
        bucketPos[0] = 0;

        for (uint32_t i = 0; i < NUM_BUCKETS - 1; i++)
        {
          bucketPos[i + 1] = bucketPos[i] + bucketSize[i];
          bucketSize[i] = 0;
        }

        uint32_t next = j + DIGIT_BITS;

        for (temp = first; temp != end; ++temp)
        {
          const size_t digit = (*temp >> j) & BUCKET_MASK;
          *(pBuffer + bucketPos[digit]++) = *temp;
          ++bucketSize[(*temp >> next) & BUCKET_MASK];
        }
      }

      end = pBuffer + (end - first);
      temp = first;
      first = pBuffer;
      pBuffer = temp;
    }
  }

  // Copy buffer to original array.
  if (first != begin)
    for (temp = first; temp != end; ++temp)
      *pBuffer++ = *temp;

  delete[] buffer;
}
