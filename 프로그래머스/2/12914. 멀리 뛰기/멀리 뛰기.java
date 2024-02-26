class Solution {
    public long solution(int n) {
        long answer = 0;
        int f0 = 0;
        int f1 = 1;
        int f2 = f0 + f1;
        for(int i = 2; i <= n; i++){
            f0 = f1;
            f1 = f2;
            f2 = (f0 + f1) % 1234567;
        }
        return f2;
    }
}