#include <stdio.h>
#include <omp.h>

long long fib_omp(int n)
{
    long long x, y;

    if (n < 2) return n;

    // n이 충분히 작을 때는 재귀 분할을 멈춰서 오버헤드 방지
    if (n < 3)
        return fib_omp(n - 1) + fib_omp(n - 2);

    // 태스크 생성
    #pragma omp task shared(x)
    x = fib_omp(n - 1);

    #pragma omp task shared(y)
    y = fib_omp(n - 2);

    // 두 태스크가 끝날 때까지 기다림
    #pragma omp taskwait
    return x + y;
}

int main(void)
{
    int n = 8;    // 계산할 피보나치 수
    long long result_omp;

    double ta_omp,tb_omp;

    omp_set_num_threads(4);

    ta_omp = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        {
            result_omp = fib_omp(n);
        }
    }
    tb_omp = omp_get_wtime();

    printf("F(%d) = %lld\n", n, result_omp);
    printf("Elapsed time = %.4f sec\n", tb_omp - ta_omp);

    return 0;
}

