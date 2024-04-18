
import java.util.*;
class Solution {
    public int[] solution(String[] maps) {
        char[][] chars = new char[maps.length][maps[0].length()];
        boolean[][] visited = new boolean[maps.length][maps[0].length()];
        List<Integer> list = new ArrayList<>();
        for(int i = 0; i < maps.length; i++) chars[i] = maps[i].toCharArray();
        for(int i = 0; i < chars.length; i++){
            for(int j = 0; j < chars[i].length; j++){
                Queue<int[]> queue = new LinkedList<>();
                queue.add(new int[] { i, j });
                int sum = 0;
                while(!queue.isEmpty()){
                    int[] q = queue.remove();
                    if(visited[q[0]][q[1]] || chars[q[0]][q[1]] == 'X') continue;
                    visited[q[0]][q[1]] = true;
                    sum += chars[q[0]][q[1]] - '0';
                    if(q[0] + 1 < chars.length) queue.add(new int[] { q[0] + 1, q[1] });
                    if(0 <= q[0] - 1) queue.add(new int[] { q[0] - 1, q[1] });
                    if(q[1] + 1 < chars[0].length) queue.add(new int[] { q[0], q[1] + 1 });
                    if(0 <= q[1] - 1) queue.add(new int[] { q[0], q[1] - 1 });
                }
                if(sum > 0) list.add(sum);
            }
        }
        if(list.size() > 0){
            Collections.sort(list);
            int[] answer = new int[list.size()];
            for(int i = 0; i < answer.length; i++) answer[i] = list.get(i);
            return answer;
        }
        return new int[] { -1 };
    }
}