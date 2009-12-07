/* 
 * prime_para.c
 *
 *
 */

#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <pthread.h>
#include <stdlib.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

typedef struct prime_arg {
  int start;
  int end;
  int sum;
  pthread_t tid;
}* prime_arg_t;

void* prime_th(void* _arg) {
  prime_arg_t arg = _arg;
  arg->sum = check_primes(arg->start, arg->end);
  return NULL;
}

int parallel_prime(int start, int end, int nthreads) {
  int t;
  prime_arg_t args = (prime_arg_t) malloc(sizeof(struct prime_arg) * nthreads);
  assert(args);
  int step = (end - start) / nthreads;
  if ((end - start) - nthreads * step != 0) step++;
  for (t = 0; t < nthreads; t++) {
    args[t].start = start + t * step;
    args[t].end = MIN(end, start + (t+1) * step);
    args[t].sum = 0;
    pthread_create(&args[t].tid, NULL, prime_th, &args[t]);
#if 0
    printf("Thread %2d: %8d - %8d : %8d\n", t, args[t].start, args[t].end, args[t].end - args[t].start);
#endif
  }
  int s = 0;
  for (t = 0; t < nthreads; t++) {
    pthread_join(args[t].tid, NULL);
    s += args[t].sum;
  }
  return s;
}

int check_prime(int n)
{
  int d;
  for (d = 2; d * d <= n; d++) {
    if (n % d == 0) return 0;
  }
#if 0
  printf("%d is prime\n", n);
#endif
  return 1;
}

int check_primes(int a, int b)
{
  int n;
  int s = 0;
  for (n = a; n < b; n++) {
    s += check_prime(n);
  }
  return s;
}

double cur_time()
{
  struct timeval tp[1];
  gettimeofday(tp, NULL);
  return tp->tv_usec * 1.0E-6 + tp->tv_sec;
}

int main(int argc, char** argv) {
  int N = 10000000;
  double t0 = cur_time();
  /*int np = check_primes(2, N); */
  if (argc < 1) {
    fprintf(stderr, "usage: %s NTHREADS", argv[0]);
    exit(1);
  }
  int npthread = atoi(argv[1]);
  int np = parallel_prime(2, N, npthread);
  double t1 = cur_time();
  printf("%d primes in %f sec\n", np, t1 - t0);
  return 0;
}

