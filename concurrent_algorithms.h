#include <algorithm> // sort, make_heap, sort_heap, swap, distance, generate, inplace_merge, minmax
#include <iostream>  // cout
#include <cstring>   // memset
#include <vector>    // used by counting sort
#include <exception> // runtime exception

using namespace concurrency;

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

#if (_MSC_VER && __cpp_lib_parallel_algorithm < 201603)
/*************************************************************************
 Concurrent MSVC++ Parallel Patterns Library radix sort.
 *************************************************************************/
template <typename T>
void pRadix(T* begin, T* end) { Concurrency::parallel_radixsort(begin, end); }

/*************************************************************************
 * Implements a concurrent MSVC++ Parallel Patterns Library quicksort.
 *************************************************************************/
template <typename T>
void ParallelQuickSortWithSTL(T* begin, T* end, int depthRemaining, long threshold)
{
  if (std::distance(begin, end) <= threshold)
    std::sort(begin, end);

  else
  {
    T* pivot = std::partition(begin, end, std::bind2nd(std::less<T>(), *begin));

    if (depthRemaining > 0)
    {
      Concurrency::parallel_invoke(
        [begin, end, pivot, depthRemaining, threshold] {
        if (pivot != end)
          ParallelQuickSortWithSTL(begin, pivot, depthRemaining - 1, threshold);
      },
        [&pivot, begin, end, depthRemaining, threshold] {
        if (pivot == begin)
          ++pivot;
        ParallelQuickSortWithSTL(pivot, end, depthRemaining - 1, threshold);
      });
    }

    else
    {
      T* pivot = std::partition(begin, end, std::bind2nd(std::less<T>(), *begin));

      if (pivot != end)
        std::sort(begin, pivot);

      if (pivot == begin)
        ++pivot;

      std::sort(pivot, end);
    }
  }
}
/*
template <typename T>
void SequentialQuickSort(T* begin, T* end, long threshold)
{
  if (std::distance(begin, end) <= threshold)
    insertionSort(begin, end);

  else
  {
    T* pivot = std::partition(begin + 1, end, std::bind2nd(std::less<int>(), *begin));

    //std::iter_swap(begin, pivot - 1);
    std::swap(*begin, *(pivot - 1));

    SequentialQuickSort(begin, pivot - 1, threshold);
    SequentialQuickSort(pivot, end, threshold);
  }
}

template <typename T>
void ParallelQuickSort(T* begin, T* end, long threshold, int depthRemaining)
{
  if (std::distance(begin, end) <= threshold)
    insertionSort(begin, end);

  else
  {
    T* pivot = std::partition(begin + 1, end, std::bind2nd(std::less<int>(), *begin));

    std::swap(*begin, *(pivot - 1));

    if (depthRemaining > 0)
    {
      Concurrency::parallel_invoke([begin, end, pivot, depthRemaining, threshold]
      {
        ParallelQuickSort(begin, pivot - 1, depthRemaining - 1, threshold);
      },
        [&pivot, begin, end, depthRemaining, threshold] {
        ParallelQuickSort(pivot, end, depthRemaining - 1, threshold);
      });
    }

    else
    {
      SequentialQuickSort(begin, pivot - 1, threshold);
      SequentialQuickSort(pivot, end, threshold);
    }
  }
}
*/

template <typename T>
void pSort(T* begin, T* end)
{
  const int maxTasks = Concurrency::CurrentScheduler::Get()->GetNumberOfVirtualProcessors();
  //ParallelQuickSort(begin, end, (int)LogN(float(maxTasks), 2.0f) + 4, 256);
  ParallelQuickSortWithSTL(begin, end, (int)LogN((float)maxTasks) + 4, 256);
}

// To generate logarithms for base 2. Log base b of a == natural log (a) / natural log (b).
inline float LogN(float value) { return logf(value) / logf(2.0f); }
#endif

#if (__cpp_lib_parallel_algorithm >= 201603)
/*************************************************************************
 * C++ STL <algorithm> library sort function with C++17 parallel option.
 * Defaults to sequential version if C++17 features not available.
 *************************************************************************/
template <typename T>
// Use parallel algorithm.
void pSort(T* begin, T* end) { std::sort(std::execution::par_unseq, begin, end); }
#endif

/*************************************************************************
 * C++ AMP radix sort function. If not available, defaults to std::sort.
*************************************************************************/
void arr_fill(concurrency::array_view<unsigned>& dest, concurrency::array_view<unsigned>& src, unsigned val)
{
  parallel_for_each(dest.extent, [dest, val, src](concurrency::index<1> idx) restrict(amp)
  {
    dest[idx] = ((unsigned)idx[0] < src.get_extent().size()) ? src[idx] : val;
  });
}

unsigned get_bits(unsigned x, unsigned numbits, unsigned bitoffset) restrict(amp) { return (x >> bitoffset) & ~(~0 << numbits); }
unsigned pow2(unsigned x) restrict(amp, cpu) { return (((unsigned)1) << x); }

unsigned tile_sum(unsigned x, concurrency::tiled_index<256> tidx) restrict(amp)
{
  unsigned l_id = tidx.local[0];
  tile_static unsigned l_sums[256][2];

  l_sums[l_id][0] = x;
  tidx.barrier.wait();

  for (unsigned i = 0; i < 8; i++)
  {
    if (l_id < pow2(7 - i))
    {
      unsigned w = (i + 1) % 2;
      unsigned r = i % 2;

      l_sums[l_id][w] = l_sums[l_id * 2][r] + l_sums[l_id * 2 + 1][r];
    }
    tidx.barrier.wait();
  }
  return l_sums[0][0];

}

unsigned tile_prefix_sum(unsigned x, concurrency::tiled_index<256> tidx, unsigned& last_val) restrict(amp)
{
  unsigned l_id = tidx.local[0];
  tile_static unsigned l_prefix_sums[256][2];

  l_prefix_sums[l_id][0] = x;
  tidx.barrier.wait();

  for (unsigned i = 0; i < 8; i++)
  {
    unsigned pow2i = pow2(i);

    unsigned w = (i + 1) % 2;
    unsigned r = i % 2;

    l_prefix_sums[l_id][w] = (l_id >= pow2i) ? (l_prefix_sums[l_id][r] + l_prefix_sums[l_id - pow2i][r]) : l_prefix_sums[l_id][r];

    tidx.barrier.wait();
  }
  last_val = l_prefix_sums[255][0];

  unsigned retval = (l_id == 0) ? 0 : l_prefix_sums[l_id - 1][0];
  return retval;
}

unsigned tile_prefix_sum(unsigned x, concurrency::tiled_index<256> tidx) restrict(amp)
{
  unsigned ll = 0;
  return tile_prefix_sum(x, tidx, ll);
}

void calc_interm_sums(unsigned bitoffset, concurrency::array<unsigned>& interm_arr, concurrency::array<unsigned>& interm_sums, concurrency::array<unsigned>& interm_prefix_sums, unsigned num_tiles)
{
  auto ext = extent<1>(num_tiles * 256).tile<256>();

  parallel_for_each(ext, [=, &interm_sums, &interm_arr](tiled_index<256> tidx) restrict(amp)
  {
    unsigned inbound = ((unsigned)tidx.global[0] < interm_arr.get_extent().size());
    unsigned num = (inbound) ? get_bits(interm_arr[tidx.global[0]], 2, bitoffset) : get_bits(0xffffffff, 2, bitoffset);

    for (unsigned i = 0; i < 4; i++)
    {
      unsigned to_sum = (num == i);
      unsigned sum = tile_sum(to_sum, tidx);

      if (tidx.local[0] == 0)
        interm_sums[i*num_tiles + tidx.tile[0]] = sum;
    }
  });

  unsigned numiter = (num_tiles / 64) + ((num_tiles % 64 == 0) ? 0 : 1);
  ext = extent<1>(256).tile<256>();

  parallel_for_each(ext, [=, &interm_prefix_sums, &interm_sums](tiled_index<256> tidx) restrict(amp)
  {
    unsigned last_val0 = 0;
    unsigned last_val1 = 0;

    for (unsigned i = 0; i < numiter; i++)
    {
      unsigned g_id = tidx.local[0] + i * 256;
      unsigned num = (g_id < (num_tiles * 4)) ? interm_sums[g_id] : 0;
      unsigned scan = tile_prefix_sum(num, tidx, last_val0);
      if (g_id < (num_tiles * 4)) interm_prefix_sums[g_id] = scan + last_val1;

      last_val1 += last_val0;
    }
  });
}

void sort_step(unsigned bitoffset, concurrency::array<unsigned>& src, concurrency::array<unsigned>& dest, concurrency::array<unsigned>& interm_prefix_sums, unsigned num_tiles)
{
  auto ext = extent<1>(num_tiles * 256).tile<256>();

  parallel_for_each(ext, [=, &interm_prefix_sums, &src, &dest](tiled_index<256> tidx) restrict(amp)
  {
    unsigned inbounds = ((unsigned)tidx.global[0] < src.get_extent().size());
    unsigned element = (inbounds) ? src[tidx.global[0]] : 0xffffffff;
    unsigned num = get_bits(element, 2, bitoffset);
    for (unsigned i = 0; i < 4; i++)
    {
      unsigned scan = tile_prefix_sum((num == i), tidx) + interm_prefix_sums[i*num_tiles + tidx.tile[0]];
      if (num == i && inbounds) dest[scan] = element;
    }
  });
}

void amp_radix_sort(concurrency::array<unsigned>& arr)
{
  unsigned size = arr.get_extent().size();
  const unsigned num_tiles = (size / 256) + ((size % 256 == 0) ? 0 : 1);
  array<unsigned> interm_arr(size);
  array<unsigned> interm_sums(num_tiles * 4);
  array<unsigned> interm_prefix_sums(num_tiles * 4);

  for (unsigned i = 0; i < 16; i++)
  {
    array<unsigned>& src = (i % 2 == 0) ? arr : interm_arr;
    array<unsigned>& dest = (i % 2 == 0) ? interm_arr : arr;
    unsigned bitoffset = i * 2;

    calc_interm_sums(bitoffset, src, interm_sums, interm_prefix_sums, num_tiles);
    sort_step(bitoffset, src, dest, interm_prefix_sums, num_tiles);
  }
}

template <typename T>
void aRadix(T* begin, T* end)
{
  accelerator default_device;
  if (default_device == accelerator(accelerator::direct3d_ref))
  {
    std::sort(begin, end);
    return;
  }
  {
    unsigned size = std::distance(begin, end);
    array<unsigned, 1> av(size, begin);

    amp_radix_sort(av);
    av.accelerator_view.wait();
    // Synchronise.
    copy(av, begin);
    av.accelerator_view.wait();
  }
}

#endif
