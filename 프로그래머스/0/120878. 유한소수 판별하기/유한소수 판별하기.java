class Solution {
    public int solution(int a, int b) {
        int d = b / gcd(a, b);
        for(int i = 2; i <= d; i++){
            boolean prime = true;
            for(int j = 2; j <= (int)Math.sqrt(i); j++){
                if(i % j == 0) prime = false;
            }
            if(prime && d % i == 0 && i != 2 && i != 5) return 2;
        }
        return 1;
    }
    
    public int gcd(int a, int b){
        if(b == 0) return a;
        return gcd(b, a % b);
    }
}