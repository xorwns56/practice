class Solution {
    public int[][] solution(int n) {
        int[][] answer = new int[n][n];
        int[] pos = new int[] {0, 0};
        int[][] direction = new int[][] {{0, 1}, {1, 0}, {0, -1}, {-1, 0}};
        int idx = 0;
        for(int i = 1; i <= n * n; i++){
            answer[pos[0]][pos[1]] = i;
            pos[0] += direction[idx][0];
            pos[1] += direction[idx][1];
            if(n <= Math.max(pos[0], pos[1]) || Math.min(pos[0], pos[1]) < 0 || answer[pos[0]][pos[1]] != 0){
                int new_idx = (idx + 1) % direction.length;
                pos[0] += direction[new_idx][0] - direction[idx][0];
                pos[1] += direction[new_idx][1] - direction[idx][1];
                idx = new_idx;
            }
        }
        return answer;
    }
}