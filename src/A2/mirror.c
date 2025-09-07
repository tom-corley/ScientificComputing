#include <stdio.h>
#include <stdlib.h>

void print_array(double *arr, int N) {
    for (int i = 0; i < N; i++) {
        printf("%.2f ", arr[i]);
    }
    printf("\n");
}

long* mirror(long *arr, long n) {
    long *ans = malloc(2 * n * sizeof(long));
    if (ans == NULL) return NULL;  

    for (long i = 0; i < n; i++) {
        ans[i] = arr[i];
        ans[2 * n - 1 - i] = arr[i];  // Set the mirrored element
    }

    return ans;
}

long int *move_sum(long int *arr, long int N) {
    long int *ans = malloc(N*sizeof(long int));
    if (ans == NULL) return NULL;
    for (long int i = 0; i < N; i++) {
        ans[i] = arr[i];
        if (i > 0) {ans[i] += arr[i - 1];}
        if (i < N - 1) {ans[i] += arr[i + 1];}
    }
    return ans;
}

long int n_elements_in_range(long int *arr, long int N, long int a, long int b) {
    long int ct = 0;
    for (int i = 0; i < N; i++) {
        if (arr[i] <= b && arr[i] >= a) {ct += 1;}
    }
    return ct;
}

int *move_to_front(int *arr, int arr_length) {
    int *ans = malloc((arr_length)*sizeof(int));
    ans[0] = arr[arr_length-1];
    for (int i = 0; i < arr_length - 1; i++) {
        ans[i] = arr[i+1];
    }
    return ans;
}

double *move_to_position(double *farr, int *iarr, int N) {
    double *ans = malloc(N*sizeof(double));
    if (ans == NULL) return NULL;
    for (int i = 0; i < N; i++) {
        ans[iarr[i]] = farr[i];
    }
    return ans;
}

int main() {
    // Test Case 4: Single Element
    double farr4[] = {3.14};
    int iarr4[] = {0};
    int N4 = 1;
    double *result4 = move_to_position(farr4, iarr4, N4);
    printf("Test Case 4 - Expected Output: 3.14\nActual Output: ");
    print_array(result4, N4);
    free(result4);

    // Test Case 5: Permutation with All Elements Identical
    double farr5[] = {7.7, 7.7, 7.7, 7.7};
    int iarr5[] = {1, 3, 0, 2};
    int N5 = 4;
    double *result5 = move_to_position(farr5, iarr5, N5);
    printf("Test Case 5 - Expected Output: 7.7 7.7 7.7 7.7\nActual Output: ");
    print_array(result5, N5);
    free(result5);

    // Test Case 6: Larger Array
    double farr6[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    int iarr6[] = {5, 3, 0, 4, 1, 2};
    int N6 = 6;
    double *result6 = move_to_position(farr6, iarr6, N6);
    printf("Test Case 6 - Expected Output: 0.3 0.5 0.6 0.2 0.4 0.1\nActual Output: ");
    print_array(result6, N6);
    free(result6);

    // Test Case 7: Duplicate Indices (testing robustness)
    double farr7[] = {4.4, 5.5, 6.6};
    int iarr7[] = {1, 1, 2};  // Duplicate index
    int N7 = 3;
    double *result7 = move_to_position(farr7, iarr7, N7);
    printf("Test Case 7 - Expected Output: Undefined behavior due to duplicate indices\nActual Output: ");
    print_array(result7, N7);
    free(result7);

    // Test Case 8: Random Permutation
    double farr8[] = {1.5, 2.5, 3.5, 4.5, 5.5};
    int iarr8[] = {4, 0, 3, 2, 1};
    int N8 = 5;
    double *result8 = move_to_position(farr8, iarr8, N8);
    printf("Test Case 8 - Expected Output: 2.5 5.5 4.5 3.5 1.5\nActual Output: ");
    print_array(result8, N8);
    free(result8);

    return 0;
}