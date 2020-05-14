#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

template<class T>
void merge(std::vector<T>& vec, int begin, int mid, int end) {
  std::vector<T> tmp(end-begin+1);
  int left = begin;
  int right = mid+1;
#pragma omp parallel for
  for (int i=0; i<tmp.size(); i++) { 
    if (left > mid)
      tmp[i] = vec[right++];
    else if (right > end)
      tmp[i] = vec[left++];
    else if (vec[left] <= vec[right])
      tmp[i] = vec[left++];
    else
      tmp[i] = vec[right++]; 
  }
  for (int i=0; i<tmp.size(); i++) 
    vec[begin++] = tmp[i];
}

template<class T>
void merge_sort(std::vector<T>& vec, int begin, int end) {
	int i, j, t;
	int length = end - begin+1;
#pragma omp parallel for private(i, t) shared(length, vec)
	//merge single ones
	for (i = 0; i < length / 2; i++)
		if (vec[i * 2] > vec[i * 2 + 1]) {
			t = vec[i * 2];
			vec[i * 2] = vec[i * 2 + 1];
			vec[i * 2 + 1] = t;
		}

	//i represents  length of merge
	for (i = 2; i < end; i *= 2) {
#pragma omp parallel for private(j) shared(end, i)
		for (j = 0; j <= end - i; j += i * 2) {
			merge(vec,j, j + i-1, (j + i * 2 -1< end ? j + i * 2-1 : end));
		}
	}
}

int main() {
  int n = 20;
  std::vector<int> vec(n);
  for (int i=0; i<n; i++) {
    vec[i] = rand() % (10 * n);
    printf("%d ",vec[i]);
  }
  printf("\n");

  merge_sort(vec, 0, n-1);

  for (int i=0; i<n; i++) {
    printf("%d ",vec[i]);
  }
  printf("\n");
}
