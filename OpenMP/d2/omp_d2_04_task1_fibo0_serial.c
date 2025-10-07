#include <stdio.h>

long long fib_serial(int n)
{
    if (n < 2) return n;
    return fib_serial(n - 1) + fib_serial(n - 2);
}

int main(void)
{
    int n = 8;    // 계산할 피보나치 수
    long long result;

    result = fib_serial(n);

    printf("F(%d) = %lld\n", n, result);

    return 0;
}


