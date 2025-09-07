#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int order_of_magnitude(double a) {
    int res = 0;
    while (a >= 10.0) {
        a /= 10.0;
        res += 1;
    }
    while (a > 0 && a < 1.0) {
        a *= 10.0;
        res -= 1;
    }

    return res;
}

int main() {
    printf("%d", order_of_magnitude(.005564));
}