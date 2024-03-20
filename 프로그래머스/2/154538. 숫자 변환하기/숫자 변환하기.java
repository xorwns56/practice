import java.util.*;
class Solution {
    public int solution(int x, int y, int n) {
        Queue<int[]> queue = new LinkedList<>();
        boolean[] visited = new boolean[y];
        queue.add(new int[] { x, 0 });
        while(!queue.isEmpty()){
            int[] tmp = queue.remove();
            if(tmp[0] == y) return tmp[1];
            else if(tmp[0] < y){
                if(!visited[tmp[0]]) {
                    visited[tmp[0]] = true;
                    queue.add(new int[] { tmp[0] + n, tmp[1] + 1 });
                    queue.add(new int[] { tmp[0] * 2, tmp[1] + 1 });
                    queue.add(new int[] { tmp[0] * 3, tmp[1] + 1 });
                }
            }
        }
        return -1;
    }
}