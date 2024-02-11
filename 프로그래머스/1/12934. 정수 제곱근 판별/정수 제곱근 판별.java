class Solution {
    public long solution(long n) {
        long x = 1;
        while(x * x < n) x++;
        return x * x == n ? (x + 1) * (x + 1) : -1;
    }
}