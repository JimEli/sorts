/*
// cocktail sort pseudocode
function cocktailSort( A : list of sortable items )
  do {
   swapped := false
   for each i in 0 to length( A ) - 2 do
   if A[ i ] > A[ i+1 ] then // test whether the two
                 // elements are in the wrong
                 // order
     swap( A[ i ], A[ i+1 ] ) // let the two elements
                // change places
     swapped := true;
   if swapped = false then
   // we can exit the outer loop here if no swaps occurred.
   break do-while loop;
   swapped := false
   for each i in length( A ) - 2 down to 0 do
   if A[ i ] > A[ i+1 ] then
     swap( A[ i ], A[ i+1 ] )
     swapped := true;
  } while swapped; // if no elements have been swapped,
                   // then the list is sorted
*/

#include <algorithm>
#include <iostream>
#include <iterator>

template <typename RandomAccessIterator>
void cocktail_sort(RandomAccessIterator begin, RandomAccessIterator end) 
{
  bool swapped = true;
  while (begin != end-- && swapped) 
  {
    swapped = false;
    for (auto i = begin; i != end; ++i) 
    {
      if (*(i + 1) < *i) 
      {
        std::iter_swap(i, i + 1);
        swapped = true;
      }
    }
    if (!swapped) 
      break;
    swapped = false;
    for (auto i = end - 1; i != begin; --i) 
    {
      if (*i < *(i - 1)) 
      {
        std::iter_swap(i, i - 1);
        swapped = true;
      }
    }
    ++begin;
  }
}
