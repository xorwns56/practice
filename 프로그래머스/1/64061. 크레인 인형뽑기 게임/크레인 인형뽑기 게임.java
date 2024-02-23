import java.util.*;
class Solution {
    public int solution(int[][] board, int[] moves) {
        int answer = 0;
        Stack<Integer> stack = new Stack<>();
        for(int i = 0; i < moves.length; i++){
            int x = moves[i] - 1;
            int y = 0;
            while(y < board.length){
                if(board[y][x] != 0){
                    if(!stack.isEmpty() && board[y][x] == (int)stack.peek()){
                        stack.pop();
                        answer += 2;
                    }else stack.push(board[y][x]);
                    board[y][x] = 0;
                    break;
                }
                y++;
            }
        }
        return answer;
    }
}