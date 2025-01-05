class Solution {
    public int solution(int n, int m, int[][] edge_list, int k, int[] gps_log) {
        int answer = 0;
        boolean[][] edge_exist = new boolean[n + 1][n + 1];
        for(int i = 1; i <= n; i++) edge_exist[i][i] = true;
        for(int i = 0; i < edge_list.length; i++) edge_exist[edge_list[i][0]][edge_list[i][1]] = edge_exist[edge_list[i][1]][edge_list[i][0]] = true;
        int[][] error = new int[k][n + 1];
        int INF = Integer.MAX_VALUE;
        for(int i = 0; i < k; i++){
            for(int j = 1; j <= n; j++){
                error[i][j] = INF;
            }
        }
        error[0][gps_log[0]] = 0;
        for(int i = 1; i < k; i++){
            int curr = gps_log[i];
            for(int j = 1; j <= n; j++){
                if(error[i - 1][j] == INF) continue;
                for(int l = 1; l <= n; l++){
                    if(!edge_exist[j][l]) continue;
                    if(l == curr) error[i][l] = Math.min(error[i][l], error[i - 1][j]);
                    else error[i][l] = Math.min(error[i][l], error[i - 1][j] + 1);;
                }
            }
        }
        return error[k - 1][gps_log[k - 1]] == INF ? -1 : error[k - 1][gps_log[k - 1]];
    }
}