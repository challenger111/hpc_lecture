#include <cstdio>
#include <cstdlib>
#include <vector>


__global__ void initial_0(int *bucket) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	bucket[i] = 0;
}

__global__ void counter(int *bucket,int *key) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	atomicAdd(&bucket[key[i], 1);
}

__global__ void recount(int *bucket, int *key) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= 50)
		return;
	for (int j = 0, k = 0; j <= i; k++) {
		key[i] = k;
		j += bucket[k];
	}
}

int main() {
	int n = 50;
	int range = 5;
	int *key;
	int *bucket;
	cudaMalloc(&key, n * sizeof(int));
	cudaMalloc(&bucket, range * sizeof(int));
	for (int i = 0; i < n; i++) {
		key[i] = rand() % range;
		printf("%d ", key[i]);
	}
	printf("\n");

	initial_0 <<<1, range >>> (bucket);
	counter <<<1, n >>> (bucket, key);
	recounter <<<(n + 32 - 1) / 32, 32, range >>> (bucket, key);


	for (int i = 0; i < n; i++) {
		printf("%d ", key[i]);
	}
	printf("\n");
}
