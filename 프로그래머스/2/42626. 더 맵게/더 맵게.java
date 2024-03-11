import java.util.*;
class Solution {
    public int solution(int[] scoville, int K) {
        PriorityQueue<Integer> queue = new PriorityQueue<>();
        int mix_count = 0;
        for(int i = 0; i < scoville.length; i++) queue.add(scoville[i]);
        while(queue.size() > 1){
            int food1 = queue.remove();
            if(food1 >= K) break;
            int food2 = queue.remove();
            queue.add(food1 + food2 * 2);
            mix_count++;
        }
        return queue.peek() < K ? -1 : mix_count;
    }
}