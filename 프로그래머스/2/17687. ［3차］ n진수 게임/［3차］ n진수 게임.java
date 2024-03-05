class Solution {
    public String solution(int n, int t, int m, int p) {
        int count = 0;
        String answer = "";
        while(count < t){
            answer += getChar(n, m * count + p - 1);
            count++;
        }
        return answer;
    }
    
    public char getChar(int n, int p){
        int digit = 1;
        int start = 0;
        int tmp;
        while(p - (tmp = ((int)Math.pow(n, digit) - start) * digit) >= 0){
            p -= tmp;
            digit++;
            start = (int)Math.pow(n, digit - 1);
        }
        int num = start + p / digit;
        char[] chars = new char[digit];
        int idx = digit - 1;
        while(idx >= 0){
            chars[idx--] = (char)(num % n + (num % n > 9 ? 'A' - 10 : '0'));
            num /= n;
        }
        return chars[p % digit];
    }
    
}