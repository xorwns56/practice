import java.util.*;
class Solution {
    public int solution(int[][] maps) {
        Queue<int[]> q = new LinkedList<int[]>();
        int[] dy = new int[] { 1, -1, 0, 0};
        int[] dx = new int[] { 0, 0, 1, -1};
        q.add(new int[]{ 0, 0 });
        while(q.size() > 0){
            int[] pos = q.remove();
            if(pos[0] + 1 == maps.length && pos[1] + 1 == maps[0].length) return maps[pos[0]][pos[1]];
            for(int i = 0; i < 4; i++){
                int[] new_pos = new int[] { pos[0] + dx[i], pos[1] + dy[i] };
                if(!(0 <= new_pos[0] && new_pos[0] < maps.length && 0 <= new_pos[1] && new_pos[1] < maps[0].length)) continue;
                else if(maps[new_pos[0]][new_pos[1]] == 0) continue;
                else if(maps[new_pos[0]][new_pos[1]] == 1){
                    maps[new_pos[0]][new_pos[1]] = maps[pos[0]][pos[1]] + 1;
                    q.add(new_pos);
                }
            }
        }
        return -1;
    }
}