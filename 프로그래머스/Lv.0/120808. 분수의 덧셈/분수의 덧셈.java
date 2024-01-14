class Solution {
    public int[] solution(int numer1, int denom1, int numer2, int denom2) {
        int denom = lcm(denom1, denom2);
        int numer = numer1 * (denom / denom1) + numer2 * (denom / denom2);
        int more = gcd(denom, numer);
        return new int[]{numer / more, denom / more};
    }
    
    public int gcd(int a, int b){
        if(b == 0) return a;
        return gcd(b, a % b);
    }
    
    public int lcm(int a, int b){
        return a * b / gcd(a, b);
    }
}