class Solution {
    public int solution(int[][] board) {
        int answer = 0;
        for(int i = 0; i < board.length; i++){
            for(int j = 0; j < board[i].length; j++){
                if(!danger(board, i, j)) answer++;
            }
        }
        return answer;
    }
    
    public boolean mine(int[][] board, int x, int y){
        return 0 <= x && x < board.length && 0 <= y && y < board[x].length && board[x][y] == 1;
    }
    
    public boolean danger(int[][] board, int x, int y){
        for(int i = -1; i <= 1; i++){
            for(int j = -1; j <= 1; j++){
                if(mine(board, x + i, y + j)) return true;
            }
        }
        return false;
    }
}