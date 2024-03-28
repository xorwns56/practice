import java.util.*;
class Solution {
    public long solution(int[] weights) {
        long answer = 0;
        Arrays.sort(weights);
        HashMap<Integer, Long> map = new HashMap<>();
        for(int i = 0; i < weights.length; i++){
            if(map.containsKey(weights[i])) answer += map.get(weights[i]);
            map.put(weights[i], map.getOrDefault(weights[i], 0L) + 1);
            if(weights[i] % 2 == 0) map.put(weights[i] / 2 * 3, map.getOrDefault(weights[i] / 2 * 3, 0L) + 1);
            if(weights[i] % 3 == 0) map.put(weights[i] / 3 * 4, map.getOrDefault(weights[i] / 3 * 4, 0L) + 1);
            map.put(weights[i] * 2, map.getOrDefault(weights[i] * 2, 0L) + 1);
        }
        return answer;
    }
}