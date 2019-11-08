template<class T>
void moveDown(T data[], int first, int last)
{
  int largest = 2 * first + 1;

  while (largest <= last)
  {
    if (largest < last && // first has 2 children (at 2*first+1 and
      data[largest] < data[largest + 1]) // 2*first+2) and the 2nd
      largest++;                         // is larger than the 1st;

    if (data[first] < data[largest]) {     // if necessary,
      swap(data[first], data[largest]);  // swap child and parent,
      first = largest;                   // and move down;
      largest = 2 * first + 1;
    }
    else
      largest = last + 1; // to exit the loop: the heap property
  }                            // isn't violated by data[first];
}

template<class T>
void heapSort(T data[], const int n)
{
  int i;

  for (i = n / 2 - 1; i >= 0; --i) // create the heap;
    moveDown(data, i, n - 1);

  for (i = n - 1; i >= 1; --i)
  {
    swap(data[0], data[i]);    // move the largest item to data[i];
    moveDown(data, 0, i - 1);  // restore the heap property;
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
 * Optimized (and efficient) heap sort function.
 * Implementation taken from:
 * http://en.wikibooks.org/wiki/Algorithm_implementation/Sorting/Heapsort
*************************************************************************/
template<typename T>
static void heapSort(T* begin, T* end)
{
  unsigned int n = std::distance(begin, end);
  unsigned int i = n / 2, parent, child;
  int t;

  // Loops until arr is sorted 
  for (;;)
  {
    // First stage - Sorting the heap
    if (i > 0)
    {
      i--;                          // Save its index to i.
      t = *(begin + i);             // Save parent value to t.
    }
    else
    {                                 // Second stage - Extracting elements in-place.
      n--;                          // Make the new heap smaller.
      if (n == 0)
        return;                   // When the heap is empty, we are done.
      t = *(begin + n);             // Save last value (it will be overwritten).
      *(begin + n) = *begin;        // Save largest value at the end of arr.
    }
    parent = i;                       // We will start pushing down t from parent.
    child = i*2 + 1;                  // Parent's left child.

    // Sift operation - pushing the value of t down the heap.
    while (child < n)
    {
      if (child + 1 < n  &&  *(begin + child + 1) > *(begin + child))
        child++;                  // Choose the largest child.

      // If any child is bigger than the parent 
      if (*(begin + child) > t)
      {
        // Move the largest child up. 
        *(begin + parent) = *(begin + child);
        parent = child;           // Move parent pointer to this child.
        child = parent*2 + 1;     // Find the next child.
      }
      else
        break;                    // t's place is found.
    }
    *(begin + parent) = t;            // We save t in the heap.
  }
}

/*************************************************************************
 * C++ STL <algorithm> library heap sort function.
*************************************************************************/
template<typename T>
void heapSort(T* begin, T* end)
{
  std::make_heap(begin, end);
  std::sort_heap(begin, end);
}
