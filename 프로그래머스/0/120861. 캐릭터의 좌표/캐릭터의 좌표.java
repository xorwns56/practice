import java.util.*;
class Solution {
    public int[] solution(String[] keyinput, int[] board) {
        HashMap<String, int[]> map = new HashMap<>();
        map.put("left", new int[]{-1, 0});
        map.put("right", new int[]{1, 0});
        map.put("up", new int[]{0, 1});
        map.put("down", new int[]{0, -1});
        int[] answer = new int[]{0, 0};
        for(int i = 0; i < keyinput.length; i++){
            int[] input = map.get(keyinput[i]);
            answer[0] = Math.max(-board[0] / 2, Math.min(board[0] / 2,  answer[0] + input[0]));
            answer[1] = Math.max(-board[1] / 2, Math.min(board[1] / 2,  answer[1] + input[1]));;
        }
        return answer;
    }
}