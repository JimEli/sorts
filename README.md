# Sorts
### Collection of sorting functions.

Generic sorting algorithms. Overloading of < and = required.
Conflict with <algorithms>, <queue>. May require a swap function: 
``` template<class T> inline void swap (T& e1, T& e2) { T tmp = e1; e1 = e2; e2 = tmp; } ```

```
template <typename T> 
void stlSort(T*, T*);       // STL sort algorithm.
```
```
template <typename T> 
void heapSort(T*, T*);      // STL heap sort algorithm.
```
```
template <typename T> 
void stlQSort(T*, T*);      // STL quicksort algorithm.
```
```
template <typename T> 
void stableSort(T*, T*);    // STL stable sort algoritm).
```
```
template <typename T> 
void insertionSort(T*, T*); // My insertion sort.
```
```
template <typename T> 
void selectionSort(T*, T*); // My selection sort algorithm.
```
```
template <typename T> 
void bubbleSort(T*, T*);
```
```
template <typename T> 
void shakerSort(T*, T*);
```
```
template <typename T> 
void combSort(T*, T*);      // My comb sort algorithm.
```
```
template <typename T> 
void myHeapSort(T*, T*);    // My heap sort algorithm.
```
```
template <typename T> 
void mergeSort(T*, T*);     // My merge sort algorithm.
```
```
template <typename T> 
void quickSort(T*, T*);     // My quicksort algorithm.
```
```
template <typename T> 
void countingSort(T*, T*);  // My counting sort algorithm.
```
```
template <typename T> 
void radixSort(T*, T*);     // My radix sort algorithm.
```

### Notes:
* (1) Insertion, selection, comb, counting, myHeap, merge and quicksort adapted from the textbook "Data Structures and Algorithms in C++", 4th edition, by Drozdek.
* (2) Radix sort inspiration from studying the Electronic Arts Standard Template Library (EASTL) source code.
* (3) Some code retrieved from here https://rosettacode.org/wiki/Category:Sorting_Algorithms
* (4) Compiled/tested with MS Visual Studio 2017 Community (v141), and Windows SDK version 10.0.17134.0
* (5) Compiled/tested with Eclipse Oxygen.3a Release (4.7.3a), using CDT 9.4.3 and MinGw gcc-g++ (6.3.0-1).
