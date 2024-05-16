class Solution {
    public int solution(int n) {
        int answer = 0;
        int count = 0;
        while(count < n){
            answer++;
            while(three(answer)) answer++;
            count++;
        }
        return answer;
    }
    public boolean three(int n){
        if(n % 3 == 0) return true;
        while(n > 0){
            if(n % 10 == 3) return true;
            n /= 10;
        }
        return false;
    }
}