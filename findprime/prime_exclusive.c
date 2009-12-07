/* 
 * prime_exclusive.c
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
  int tid;
  int stolen;
  int jobs;
}* prime_arg_t;

pthread_mutex_t m;

void* prime_th(void* _arg) {
  prime_arg_t arg = _arg;
  arg->sum = check_primes(arg->start, arg->end);
  return NULL;
}

void* do_job(void* _arg) {
  prime_arg_t arg = _arg;
  int j = 0;
#if 0
  printf("JOBS = %d\n", arg[0].jobs);
#endif
  while(j < arg[0].jobs) {
    pthread_mutex_lock(&m);
    if (arg[j].stolen == 0) {
      arg[j].stolen = 1;
      pthread_mutex_unlock(&m);
      arg[j].sum = check_primes(arg[j].start, arg[j].end);
    } else {
      pthread_mutex_unlock(&m);
    }
    j++;
  }
  return NULL;
}

parallel_prime_mutex(int start, int end, int nthreads, int divide) {
  int t,  j;
  prime_arg_t args = (prime_arg_t) malloc(sizeof(struct prime_arg) * divide);
  pthread_t* threads = (pthread_t*) malloc(sizeof(pthread_t) * nthreads);
  assert(threads);
  assert(args);
  int step = (end - start) / divide;
  if ((end - start) - divide * step != 0) step++;
  for (j = 0; j < divide; j++) {
    args[j].start = start + j * step;
    args[j].end = MIN(end, start + (j+1) * step);
    args[j].sum = 0;
    args[j].stolen = 0;
    args[j].jobs = divide;
#if 0
    printf("Job %4d: %10d - %10d : %10d\n", j, args[j].start, args[j].end, args[j].end - args[j].start);
#endif
  }
  for (t = 0; t < nthreads; t++) {
    pthread_create(&threads[t], NULL, do_job, args);
  }
  for (t = 0; t < nthreads; t++) {
    pthread_join(threads[t], NULL);
    printf("JOIN %d(%ld)\n", t, threads[t]);
  }
  int s = 0;
  for (j = 0; j < divide; j++) {
    s += args[j].sum;
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
  if (argc != 3) {
    fprintf(stderr, "usage: %s NTHREADS JOBS", argv[0]);
    exit(1);
  }
  int npthread = atoi(argv[1]);
  int jobs = atoi(argv[2]);
  int np = parallel_prime_mutex(2, N, npthread, jobs);
  double t1 = cur_time();
  printf("%d primes in %f sec\n", np, t1 - t0);
  return 0;
}

