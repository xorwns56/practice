class Solution {
    public int solution(int n) {
        return nextBigNumber(n);
    }
    public int nextBigNumber(int n) {
        int postPattern = n & -n;
        int smallPattern = ((n ^ (n + postPattern)) / postPattern) >> 2;
        return n + postPattern | smallPattern;
    }
}