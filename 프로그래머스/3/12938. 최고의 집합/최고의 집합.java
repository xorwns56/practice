class Solution {
    public int[] solution(int n, int s) {
        int[] answer = new int[n];
        for(int i = 0; i < answer.length; i++) answer[i] = s / n + (i < answer.length - s % n ? 0 : 1);
        return answer[0] == 0 ? new int[] { -1 } : answer;
    }
}