void inline clear(long& q) { q = -1; }
int inline isEmpty(long q) { return q == -1; }

template<class T>
class RadixSortNode 
{
public:
  T *arr;
  RadixSortNode *next;
};

void radixSort2(long data[], const int n) 
{
  register int d, j, k, factor, where;
  RadixSortNode<long> n1, n2, *p1, *p2;

  n1.arr = data;
  n2.arr = new long[n];

  for (j = 0; j < n; j++)
    n2.arr[j] = data[j];

  long *queue = new long[n], queueHeads[radix], queueTails[radix];

  p1 = n2.next = &n1;
  p2 = n1.next = &n2;

  for (d = 0, factor = 1; d < digits; factor *= radix, d++) 
  {
    for (j = 0; j < radix; j++)
      clear(queueHeads[j]);
    
    for (j = 0; j < n; j++) 
    {
      where = (p1->arr[j] / factor) % radix;

      if (isEmpty(queueHeads[where]))
        queueTails[where] = queueHeads[where] = j;
      else 
      {
        queue[queueTails[where]] = j;
        queueTails[where] = j;
      }
    }
    
    for (j = 0; j < radix; j++)
      if (!(isEmpty(queueHeads[j])))
        clear(queue[queueTails[j]]);
    
    for (j = k = 0; j < radix; j++)
      while (!(isEmpty(queueHeads[j]))) 
      {
        p2->arr[k++] = p1->arr[queueHeads[j]];
        queueHeads[j] = queue[queueHeads[j]];
      }

    p2 = p2->next;
    p1 = p1->next;
  }

  if (digits % 2 != 0) // if digits is an odd number;
    for (d = 0; d < n; d++)
      data[d] = p1->arr[d];
}

class RadixSortNode 
{
public:
  long *arr;
  RadixSortNode *next;
  RadixSortNode() { next = 0; }
  
  RadixSortNode(long *a, int n) 
  {
    arr = new long[n];

    for (int i = 0; i < n; i++)
      arr[i] = a[i];
    next = 0;
  }
  
  RadixSortNode(int n) 
  {
    arr = new long[n];
    next = 0;
  }
};

void bitRadixSort2(long data[], const int n, int b) 
{
  int pow2b = 1;
  pow2b <<= b;
  int d, j, k, where, pos = 0, mask = pow2b - 1;
  int last = (bits % b == 0) ? (bits / b) : (bits / b + 1);
  long *queues = new long[n], *queueHeads = new long[pow2b];
  long *queueTails = new long[pow2b];
  RadixSortNode *n2 = new RadixSortNode(data, n), *n1 = new RadixSortNode(n);

  n1->arr = data;
  n2->next = n1;
  n1->next = n2;
  for (d = 0; d < last; d++) 
  {
    for (j = 0; j < pow2b; j++)
      clear(queueHeads[j]);
    
    for (j = 0; j < n; j++) 
    {
      where = (n1->arr[j] & mask) >> pos;

      if (isEmpty(queueHeads[where]))
        queueTails[where] = queueHeads[where] = j;
      else 
      {
        queues[queueTails[where]] = j;
        queueTails[where] = j;
      }
    }
    
    mask <<= b;
    pos = pos + b;
    for (j = 0; j < pow2b; j++)
      if (!(isEmpty(queueHeads[j])))
        clear(queues[queueTails[j]]);
    
    for (j = k = 0; j < pow2b; j++)
      while (!(isEmpty(queueHeads[j]))) 
      {
        n2->arr[k++] = n1->arr[queueHeads[j]];
        queueHeads[j] = queues[queueHeads[j]];
      }

    n2 = n2->next;
    n1 = n1->next;
  }

  if (last % 2 != 0) // if bits is an odd number;
    for (d = 0; d < n; d++)
      data[d] = n1->arr[d];
}
