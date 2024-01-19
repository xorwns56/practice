class Solution {
    public int solution(int n) {
        double sqrt = Math.sqrt(n);
        return (int)sqrt == (int)Math.ceil(sqrt) ? 1 : 2;
    }
}