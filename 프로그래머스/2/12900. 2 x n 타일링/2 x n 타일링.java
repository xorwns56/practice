class Solution {
    public int solution(int n) {
        if(n <= 2) return n;
        int[] f = new int[3];
        f[0] = 1;
        f[1] = 1;
        f[2] = f[0] + f[1];
        for(int i = 3; i <= n; i++){
            f[0] = f[1];
            f[1] = f[2];
            f[2] = (f[0] + f[1]) % 1000000007;
        }
        return f[2];
    }
}