import java.util.*;
class Solution {
    public int solution(int a, int b, int c, int d) {
        PriorityQueue<int[]> queue = new PriorityQueue<>((arr1, arr2)->{ return arr1[1] != arr2[1] ? arr2[1] - arr1[1] : 1; });
        int[] dice = new int[7];
        dice[a]++; dice[b]++; dice[c]++; dice[d]++;
        for(int i = 1; i <= 6; i++){
            if(dice[i] > 0) queue.add(new int[]{i, dice[i]});
        }
        if(queue.size() == 1) return queue.remove()[0] * 1111;
        else if(queue.size() == 2){
            int[] p = queue.remove();
            int[] q = queue.remove();
            return p[1] == 3 ? (10 * p[0] + q[0]) * (10 * p[0] + q[0]) : ((p[0] + q[0]) * Math.abs(p[0] - q[0]));
        }else if(queue.size() == 3){
            queue.remove();
            return queue.remove()[0] * queue.remove()[0];
        }else if(queue.size() == 4) return Math.min(Math.min(queue.remove()[0], queue.remove()[0]), Math.min(queue.remove()[0], queue.remove()[0]));
        return 0;
    }
}