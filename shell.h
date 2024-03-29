template<class T>
void shellSort(T data[], const int n)
{
  register int i, j, hCnt, h;
  int increments[20], k;

  // create an appropriate number of increments h
  for (h = 1, i = 0; h < n; i++) 
  {
    increments[i] = h;
    h = 3 * h + 1;
  }
  
  // loop on the number of different increments h
  for (i--; i >= 0; i--) 
  {
    h = increments[i];
    // loop on the number of subarrays h-sorted in ith pass
    for (hCnt = h; hCnt < 2 * h; hCnt++) 
    {
      // insertion sort for subarray containing 
      // every hth element of array data
      for (j = hCnt; j < n;) 
      {
        T tmp = data[j];

        k = j;
        while (k - h >= 0 && tmp < data[k - h]) 
        {
          data[k] = data[k - h];
          k -= h;
        }
        data[k] = tmp;
        j += h;
      }
    }
  }
}
