class Solution {
    public int solution(int n) {
        int answer = 0;
        int i = 1;
        while(n > 0){
            if(n % i == 0) answer++;
            n -= i;
            i++;
        }
        return answer;
    }
}