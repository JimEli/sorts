/*
// pseudocode
function combsort(array input)
  gap := input.size //initialize gap size
  loop until gap = 1 and swaps = 0
    //update the gap value for a next comb. Below is an example
    gap := int(gap / 1.25)
    if gap < 1
      //minimum gap is 1
      gap := 1
    end if
    i := 0
    swaps := 0 //see Bubble Sort for an explanation
    //a single "comb" over the input list
    loop until i + gap >= input.size //see Shell sort for similar idea
      if input[i] > input[i+gap]
        swap(input[i], input[i+gap])
        swaps := 1 // Flag a swap has occurred, so the
               // list is not guaranteed sorted
      end if
      i := i + 1
    end loop
  end loop
end function
*/

template <class T>
void combSort(T data[], const int n)
{
  int step = n, j, k;
  
  while ((step = int(step / 1.3)) > 1)        // phase 1
    for (j = n - 1; j >= step; j--) 
    {
      k = j - step;
      if (data[j] < data[k])
        swap(data[j], data[k]);
    }
  
  bool again = true;
  
  for (std::size_t i = 0; i < n - 1 && again; i++)    // phase 2
    for (j = n - 1, again = false; j > i; --j)
      if (data[j] < data[j - 1]) 
      {
        swap(data[j], data[j - 1]);
        again = true;
      }
}
