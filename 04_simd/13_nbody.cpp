#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
    int jarr[N];   // count of j
    for(int k=0;k<N;k++)
      jarr[k]=k;
    __m256 xivec=_mm256_set1_ps(x[i]);
    __m256 yivec=_mm256_set1_ps(y[i]);
    __m256 xjvec=_mm256_load_ps(x);
    __m256 yjvec=_mm256_load_ps(y);
    __m256 rxvec=_mm256_sub_ps(xivec,xjvec);   //float rx = x[i] - x[j];
    __m256 ryvec=_mm256_sub_ps(yivec,yjvec);   //float ry = y[i] - y[j];
    __m256 t1vec=_mm256_mul_ps(rxvec,rxvec);
    __m256 t2vec=_mm256_mul_ps(ryvec,ryvec);
    __m256 rtvec=_mm256_add_ps(t1vec,t2vec);   //(rx * rx + ry * ry)
    
    __m256i jvec=_mm256_set1_si256(i);         
    __m256i iarrvec=_mm256_load_si256((__m256i*)iarr);
    __m256 mask=_mm256_cmpeq_ps(jvec,jarrvec);  
    __m256 rvec=_mm256_rsqrt_ps(rtvec);
    __m256 r3vec=_mm256_mul_ps(rvec,_mm256_mul_ps(rvec,rvec))
    __m256 mvec=_mm256_load_ps(m);
    //if mask element is 0,means i=j,minus value multiply to 0
    __m256 fxminus=_mm256_mul_ps(mask,_mm256_mul_ps(rxvec,_mm256_mul_ps(mvec,r3vec)));
    __m256 fyminus=_mm256_mul_ps(mask,_mm256_mul_ps(ryvec,_mm256_mul_ps(mvec,r3vec)));
    float fxm[N], fym[N];
    _mm256_store_ps(fxm, fxminus);
    _mm256_store_ps(fym, fyminus);
    
    for(int j=0;j<N;j++){
      fx[i]-=fxm[j];
      fy[i]-=fym[j];
    }
    
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
