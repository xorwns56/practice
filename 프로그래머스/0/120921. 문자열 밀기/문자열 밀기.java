class Solution {
    public int solution(String A, String B) {
        char[] chars_A = A.toCharArray();
        char[] chars_B = B.toCharArray();
        int push = 0;
        while(push < chars_A.length) {
            boolean eq = true;
            for(int i = 0; i < chars_A.length; i++){
                if(chars_A[i] != chars_B[(i + push) % chars_B.length]) eq = false;
            }
            if(eq) return push;
            push++;
        }
        return -1;
    }
}