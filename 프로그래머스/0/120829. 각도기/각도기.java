class Solution {
    public int solution(int angle) {
        int answer = 0;
        if(0 < angle) answer++;
        if(90 <= angle) answer++;
        if(90 < angle) answer++;
        if(180 <= angle) answer++;
        return answer;
    }
}