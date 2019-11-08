/*************************************************************************
* Title: Benchmark Sort Functions.
* File: sort_funcs.h
* Author: James Eli
* Date: 11/16/2018
*
* Collection of several sort functions:
*
* Notes:
*  (1) Insertion, selection, comb, counting, myHeap, merge and quicksort
*      adapted from the textbook "Data Structures and Algorithms in C++",
*      4th edition, by Drozdek.
*  (2) Radix sort inspiration from studying the Electronic Arts Standard
*      Template Library (EASTL) source code.
*  (3) Compiled/tested with MS Visual Studio 2017 Community (v141), and
*      Windows SDK version 10.0.17134.0
*  (4) Compiled/tested with Eclipse Oxygen.3a Release (4.7.3a), using
*      CDT 9.4.3 and MinGw gcc-g++ (6.3.0-1).
*
*************************************************************************
* Change Log:
*   11/16/2018: Initial release. JME
*************************************************************************/
#ifndef _SORT_FUNCTIONS_H_
#define _SORT_FUNCTIONS_H_

#include <algorithm> // sort, make_heap, sort_heap, swap, distance, generate, inplace_merge, minmax
#include <iostream>  // cout
#include <cstring>   // memset
#include <vector>    // used by counting sort
#include <exception> // runtime exception

/*
// Integral sort function prototypes.
template <typename T> void stlSort(T*, T*);       // STL sort algorithm.
template <typename T> void heapSort(T*, T*);      // STL heap sort algorithm.
template <typename T> void stlQSort(T*, T*);      // STL quicksort algorithm.
template <typename T> void stableSort(T*, T*);    // STL stable sort algoritm).
template <typename T> void insertionSort(T*, T*); // My insertion sort.
template <typename T> void selectionSort(T*, T*); // My selection sort algorithm.
template <typename T> void bubbleSort(T*, T*);
template <typename T> void shakerSort(T*, T*);
template <typename T> void combSort(T*, T*);      // My comb sort algorithm.
template <typename T> void myHeapSort(T*, T*);    // My heap sort algorithm.
template <typename T> void mergeSort(T*, T*);     // My merge sort algorithm.
template <typename T> void quickSort(T*, T*);     // My quicksort algorithm.
template <typename T> void countingSort(T*, T*);  // My counting sort algorithm.
template <typename T> void radixSort(T*, T*);     // My radix sort algorithm.
*/

/*************************************************************************
 * C++ STL <algorithm> library sort function. sort() uses an introsort
 * (quicksort which switches to heapsort when the recursion reaches a
 * certain depth).
*************************************************************************/
template<typename T>
void stlSort(T* begin, T* end) { std::sort(begin, end); }

/*************************************************************************
 * C++ STL <algorithm> library heap sort function.
*************************************************************************/
template<typename T>
void heapSort(T* begin, T* end)
{
  std::make_heap(begin, end);
  std::sort_heap(begin, end);
}

/*************************************************************************
 * C++ STL <algorithm> library qsort function.
*************************************************************************/
template<typename T>
T compare(const void *a, const void *b) { return (*(T*)a - *(T*)b); }
template<typename T>
void stlQSort(T* begin, T* end) { std::qsort(begin, std::distance(begin, end), sizeof(T), compare); }

/*************************************************************************
 * C++ STL <algorithm> library stable sort (merge) function.
*************************************************************************/
template<typename T>
void stableSort(T* begin, T* end) { std::stable_sort(begin, end); }

/*************************************************************************
 * Implements the insertion sort algorithm.
 *************************************************************************/
template<class T>
void insertionSort(T* begin, T* end)
{
  for (T* i = begin; i != end; i++)
  {
    T tmp = *i, *j;

    for (j = i; j >= begin && tmp < *(j - 1); j--)
      *j = *(j - 1);
    *j = tmp;
  }
}
/*
// Binary insertion sort O(n^2) time O(1) space, binary search to insert instead of swapping
template <typename T>
void insertion_sort(T* begin, T* end, T* cur) {
  for (; cur < end; ++cur)
  {
    T pivot = std::move(*cur);

    // upper_bound is binary search for correct insertion position
    const T pos = std::upper_bound(begin, cur, pivot);
    for (T* ins = cur; ins != pos; --ins)
      *ins = std::move(*(ins - 1));
    *pos = std::move(pivot);
  }
}
template <typename T>
void insertionSort(T* begin, T* end) { insertion_sort(begin, end, begin); }
*/

/*************************************************************************
 * Implements the selection sort algorithm.
 *************************************************************************/
template <typename T>
void selectionSort(T* begin, T* end)
{
  T *j, least;

  for (T* i = begin; i != end; i++)
  {
    for (j = (i + 1); j != end; j++)
      least = *i;

    if (*j < least)
      least = *j;

    std::swap(least, *i);
  }
}

/*************************************************************************
 * Implements the bubble sort algorithm.
 *************************************************************************/
template <typename T>
void bubbleSort(T* begin, T* end)
{
  for (T* i = begin; i < end; i++)
    for (T* j = end - 1; j != i; --j)
      if (*j < *(j - 1))
        std::swap(*j, *(j - 1));
}

/*************************************************************************
 * Implements the cocktail shaker sort algorithm.
 *************************************************************************/
template <typename T>
void shakerSort(T* begin, T* end)
{
  bool swapped = true;
  T* last = end;

  while (begin != last-- && swapped)
  {
    swapped = false;

    for (T* i = begin; i != last; ++i)
    {
      if (*(i + 1) < *i)
      {
        std::swap(*i, *(i + 1));
        swapped = true;
      }
    }

    if (!swapped)
      break;

    swapped = false;

    for (T* i = last - 1; i != begin; --i)
    {
      if (*i < *(i - 1))
      {
        std::swap(*i, *(i - 1));
        swapped = true;
      }
    }

    ++begin;
  }
}

/*************************************************************************
 * Implements a comb sort algorithm.
 *************************************************************************/
template <typename T>
void combSort(T* begin, T* end)
{
  T *current, *next;
  T length = std::distance(begin, end), space = length;

  for (bool swapped = false; (space > 1) || swapped;)
  {
    space = ((space * 10) + 3) / 13;

    if ((space == 9) || (space == 10))
      space = 11;

    current = next = begin;
    next += space;

    for (swapped = false; next != end; current++, next++)
      if (*next < *current)
      {
        std::swap(*current, *next);
        swapped = true;
      }
  }
}

/*************************************************************************
 * Implements an alternative heap sort algorithm.
 *************************************************************************/
template<class T>
void moveDown(T data[], int first, int last)
{
  int largest = 2 * first + 1;

  while (largest <= last)
  {
    if (largest < last && data[largest] < data[largest + 1])
      largest++;

    if (data[first] < data[largest])
    {
      std::swap(data[first], data[largest]);
      first = largest;
      largest = 2 * first + 1;
    }
    else
      largest = last + 1;
  }
}

template<class T>
void myHeapSort(T* begin, T* end)
{
  int n = std::distance(begin, end);

  for (int i = n / 2 - 1; i >= 0; --i)
    moveDown(begin, i, n - 1);

  for (int i = n - 1; i >= 1; --i)
  {
    std::swap(*begin, *(begin + i));
    moveDown(begin, 0, i - 1);
  }
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

/*************************************************************************
 * Implements a basic quick sort algorithm.
 *************************************************************************/
template<class T>
void quick_sort(T* begin, int first, int last)
{
  int lower = first + 1, upper = last;

  //std::swap(*(begin + first), *(begin + (first + last)/2));
  T bound = *(begin + first);

  while (lower <= upper)
  {
    while (*(begin + lower) < bound)
      lower++;

    while (bound < *(begin + upper))
      upper--;

    //if (lower < upper)
      //std::swap(*(begin + lower++), *(begin + upper--));
    //else
      //lower++;
    if (lower < upper)
    {
      std::swap(*(begin + lower), *(begin + upper));
      upper--;
    }
    lower++;
  }

  std::swap(*(begin + upper), *(begin + first));

  if (first < upper - 1)
    quick_sort(begin, first, upper - 1);

  if (upper + 1 < last)
    quick_sort(begin, upper + 1, last);
}

template<class T>
void quickSort(T* begin, T* end)
{
  const int n = std::distance(begin, end);

  if (n < 2)
    return;

  int max = 0;
  // Find largest element and put it at the end.
  for (int i = 1; i < n; i++)
    if (*(begin + max) < *(begin + i))
      max = i;

  // Largest element is now in its final position.
  std::swap(*(begin + n - 1), *(begin + max));
  quick_sort(begin, 0, n - 2);
}

/*************************************************************************
 * Implements a radix sort algorithm.
 *************************************************************************/
template <typename T>
void radixSort(T* begin, T* end)
{
  T* temp, *first = begin;
  T* buffer = new T[std::distance(begin, end)];
  T* pBuffer = buffer;
  constexpr int DIGIT_BITS = 8;                   // Bit per digit.
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

/*************************************************************************
 * Implements least significant digit radix sort algorithm.
 *************************************************************************/
 // Radix sort comparator for 32-bit two's complement integers.
class radixTest {
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

/*************************************************************************
 * Implements a parallel sort algorithm.
 *************************************************************************/
 // Count the number of elements in each bucket
 // list: array pointer
 // sz: length of array
 // np: number of processors (threads)
 // range: range of buckets
 // bucket: number of items in each bucket (return value)
template <typename T>
static void count(const T* list, size_t sz, unsigned np, const T* range, size_t* bucket)
{
  for (unsigned n = 0; n < np; n++)
    bucket[n] = 0;

  for (size_t i = 0; i < sz; i++)
  {
    unsigned n;

    for (n = 0; n < np - 1; n++)
    {
      if (list[i] < range[n])
      {
        bucket[n]++;
        break;
      }
    }

    if (n == np - 1)
      bucket[np - 1]++;
  }
}

// Reorder array to separate buckets.
// list: array pointer
// sz: length of array
// np: number of processors (threads)
// range: range of buckets
// map: reordering map
// out: reordered array (return value)
template <typename T>
static void reOrder(const T* list, size_t sz, unsigned np, const T* range, size_t* map, T* out)
{
  for (size_t i = 0; i < sz; i++)
  {
    unsigned n;

    for (n = 0; n < np - 1; n++)
    {
      if (list[i] < range[n])
      {
        out[map[n]++] = list[i];
        break;
      }
    }

    if (n == np - 1)
      out[map[np - 1]++] = list[i];
  }
}

// Sort an array using std::sort, sz: length of array, list: array pointer
template <typename T> static void _sort(T* list, size_t sz) { std::sort(list, list + sz); }

// Parallel sort
template <typename T>
void pSort(T* begin, T* end)
{
  int sz = static_cast<int>(std::distance(begin, end));
  unsigned np = 6;   // number of processors (threads).
  unsigned sf = 100; // sampling factor (100 should be good enough for most cases).
  unsigned nSample = np * sf;
  unsigned nMap = np * np;
  std::vector<T> range(np);
  std::vector<T> sample(nSample);
  std::vector<size_t> tmp(nMap, 0), map(nMap, 0), bucket(nMap, 0);
  std::vector<std::thread> threads(np);

  // Step 1: Pick np*sf random elements from the input array.
  // We use these elements to split the array into even partitions.
  for (unsigned i = 0; i < nSample; i++)
    sample[i] = begin[rand() % sz];

  std::sort(std::execution::par_unseq, &sample[0], &sample[0] + nSample);

  for (unsigned i = 0; i < np - 1; i++)
    range[i] = sample[(i + 1)*sf];

  // Step 2: Count how many elements go into each bucket. This is done in parallel.
  unsigned part = sz / np;

  for (unsigned i = 0; i < np; i++)
  {
    size_t start = i * part;
    size_t length = (i + 1 == np) ? sz - i * part : part;

    threads[i] = std::thread(count<T>, begin, length, np, &range[0], &bucket[i*np]);
  }

  for (auto&& thread : threads)
    thread.join();

  // Step 3: Create a mapping scheme to separate partitions.
  // Since the number of operations is small, we do this step in serial.
  for (unsigned i = 0; i < nMap; i++)
    tmp[i] = i ? tmp[i - 1] + bucket[((i - 1) % np)*np + (i - 1) / np] : 0;

  for (unsigned i = 0; i < nMap; i++)
    map[i] = tmp[(i%np)*np + i / np];

  for (unsigned i = 0; i < nMap; i++)
    tmp[i] = map[i];

  // Step 4: Map the input array into np separate partitions.    
  std::vector<T> sorted(sz);

  for (unsigned i = 0; i < np; i++)
  {
    size_t start = i * part;
    size_t length = (i + 1 == np) ? sz - i * part : part;

    threads[i] = std::thread(reOrder<T>, begin, length, np, &range[0], &tmp[i*np], &sorted[0]);
  }

  for (auto&& thread : threads)
    thread.join();

  // Step 5: Sort each partition in a separate thread.    
  for (unsigned i = 0; i < np; i++)
  {
    size_t start = map[i];
    size_t length = (i + 1 == np) ? sz - map[i] : map[i + 1] - map[i];

    threads[i] = std::thread(_sort<T>, &sorted[start], length);
  }

  for (auto&& thread : threads)
    thread.join();

  // Step 6: Copy the temporary container for the sorted array into the input array.    
  std::memcpy(begin, &sorted[0], sz * sizeof(T));
}

#endif

