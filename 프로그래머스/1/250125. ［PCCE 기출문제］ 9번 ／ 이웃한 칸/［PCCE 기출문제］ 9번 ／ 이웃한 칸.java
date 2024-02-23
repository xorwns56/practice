class Solution {
    public int solution(String[][] board, int h, int w) {
        int answer = 0;
        for(int i = -1; i <= 1; i++){
            for(int j = -1; j <= 1; j++){
                if(Math.abs(i) == Math.abs(j)) continue;
                int x = h + i;
                if(!(0 <= x && x < board.length)) continue;
                int y = w + j;
                if(!(0 <= y && y < board[0].length)) continue;
                answer += board[h][w].equals(board[x][y]) ? 1 : 0;
            }
        }
        return answer;
    }
}