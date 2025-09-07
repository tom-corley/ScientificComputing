int isprime(long k) {
    if (k < 2) return 0;
    for (int i = 2; i < k; i++) {
        if (k % i == 0) return 0;
    }
    return 1;
}