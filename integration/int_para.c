/* 
 * int_seq.c
 *
 * compile:
 *   gcc int_seq.c -lm
 *
 */

#include <assert.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

typedef struct integrate_arg
{
  double (*f)(double);
  double a;
  double b;
  int N;
  double s;
  pthread_t tid;
} * integrate_arg_t;

double integrate(double (*f)(double), double a, double b, int N);

void * integrate_th(void * _arg)
{
  integrate_arg_t arg = _arg;
  arg->s = integrate(arg->f, arg->a, arg->b, arg->N);
  return NULL;
}

double parallel_integrate(double (*f)(double), double a, double b, int N, 
			  int n_threads)
{
  int t;
  integrate_arg_t args = (integrate_arg_t)malloc(sizeof(struct integrate_arg) * n_threads);
  assert(args);
  for (t = 0; t < n_threads; t++) {
    args[t].f = f;
    args[t].a = a + (b - a) * t / n_threads;
    args[t].b = a + (b - a) * (t + 1) / n_threads;
    args[t].N = (N * (t + 1)) / n_threads - (N * t) / n_threads;
    pthread_create(&args[t].tid, NULL, integrate_th, &args[t]);
  }
  double s = 0.0;
  for (t = 0; t <n_threads; t++) {
    pthread_join(args[t].tid, NULL);
    s += args[t].s;
  }
  return s;
}

double integrate(double (*f)(double), double a, double b, int N)
{
  int i;
  double s = 0.0;
  for (i = 0; i < N; i++) {
    double x      = a + (b - a) *  i      / N;
    double x_next = a + (b - a) * (i + 1) / N;
    double dx = x_next - x;
    s += f(x) * dx;
  }
  return s;
}


double cur_time()
{
  struct timeval tp[1];
  gettimeofday(tp, NULL);
  return tp->tv_usec * 1.0E-6 + tp->tv_sec;
}

int main(int argc, char ** argv)
{
  double a = 0.0;
  double b = M_PI / 2.0;
  int nthreads = 2;
  if (argc > 1) nthreads = atoi(argv[1]);
  double t0 = cur_time();
  double s = parallel_integrate(sin, a, b, 100000000, nthreads);
  double t1 = cur_time();
  printf("\\int_[%f,%f] sin(x) dx = %f\n", a, b, s);
  printf("%.6f sec\n", t1 - t0);
  return 0;
}
