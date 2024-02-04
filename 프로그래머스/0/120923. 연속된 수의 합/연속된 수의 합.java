class Solution {
    public int[] solution(int num, int total) {
        int[] answer = new int[num];
        int x = (total - num * (num - 1) / 2) / num;
        for(int i = 0; i < answer.length; i++) answer[i] = x + i;
        return answer;
    }
}