class Solution {
    int answer = 0;
    public int solution(int k, int[][] dungeons) {
        dungeon(k, dungeons, new boolean[dungeons.length], 0);
        return answer;
    }
    public void dungeon(int k, int[][] dungeons, boolean[] visit, int visit_count){
        for(int i = 0; i < dungeons.length; i++){
            if(dungeons[i][0] <= k && !visit[i]){
                visit[i] = true;
                dungeon(k - dungeons[i][1], dungeons, visit, visit_count + 1);
                visit[i] = false;
            }
        }
        answer = Math.max(answer, visit_count);
    }
}