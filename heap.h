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
