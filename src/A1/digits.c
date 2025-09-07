#include <stdio.h>

int num_conseq_digits(long k) {
    int res = 1; 
    int cur_len = 1;
    int last_digit = k % 10; // Initialize with the last digit
    k /= 10;
    while (k > 0) {
        int cur_digit = k % 10;
        if (cur_digit == last_digit) {
            cur_len++;
        } else {
            if (cur_len > res) {
                res = cur_len;
            }
            cur_len = 1; // Reset the current length
        }
        last_digit = cur_digit; // Update last_digit for the next comparison
        k /= 10;
    } if (cur_len > res) {
        res = cur_len;
    }
    return res;
}

int main() {
    printf("%ld", num_conseq_digits(55555));
    return 0;
}