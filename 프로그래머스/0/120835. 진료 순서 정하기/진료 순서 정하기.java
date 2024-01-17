import java.util.*;
class Solution {
    public int[] solution(int[] emergency) {
        TreeMap<Integer, Integer> treeMap = new TreeMap<>((a, b)->a == b ? 0 : (a > b ? -1 : 1));
        for(int i = 0; i < emergency.length; i++) treeMap.put(emergency[i], i);
        int[] answer = new int[emergency.length];
        int idx = 1;
        for(Map.Entry<Integer, Integer> entry : treeMap.entrySet()) answer[entry.getValue()] = idx++;
        return answer;
    }
}