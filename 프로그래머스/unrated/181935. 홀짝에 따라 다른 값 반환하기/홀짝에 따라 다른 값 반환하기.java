class Solution {
    public int solution(int n) {
        int answer = 0;
        for(int i = n; i >= 0; i -= 2){
            answer += (i & 1) == 0 ? i * i : i;
        }
        return answer;
    }
}