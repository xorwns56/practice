class Solution {
    public int[] solution(int[][] score) {
        int[] answer = new int[score.length];
        for(int i = 0; i < answer.length; i++) answer[i] = 1;
        for(int i = 0; i < score.length; i++){
            for(int j = 0; j < score.length; j++){
                if(i == j) continue;
                if((score[i][0] + score[i][1]) / 2. > (score[j][0] + score[j][1]) / 2.) answer[j]++;
            }
        }
        return answer;
    }
}