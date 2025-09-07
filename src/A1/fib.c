#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main() {
    int res = 0;;
    int cur = 1;
    int prev = 0;
    int temp;
    int N;
    scanf("%d", &N);
    for (int i = 1; i <= N; i++) {
        res += prev;
        temp = cur;
        cur = prev + cur;
        prev = temp;
    }
    printf("%d", res);
}