class Solution {
    public int solution(int a, int b, int c) {
        int score = (a == b ? 1 : 0) + (b == c ? 1 : 0) + (a == c ? 1 : 0);
        return (a + b + c) * (1 <= score ? (a * a + b * b + c * c) : 1) * (3 <= score ? (a * a * a + b * b * b + c * c * c) : 1);
    }
}