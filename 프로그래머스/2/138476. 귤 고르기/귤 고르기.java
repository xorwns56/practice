import java.util.*;
class Solution {
    public int solution(int k, int[] tangerine) {
        HashMap<Integer, Integer> map = new HashMap<>();
        for(int i = 0; i < tangerine.length; i++) map.put(tangerine[i], map.getOrDefault(tangerine[i], 0) + 1);
        Integer[] arr = map.values().toArray(new Integer[0]);
        Arrays.sort(arr);
        int answer = 0;
        for(int i = arr.length - 1; i >= 0 && k > 0; i--){
            k -= arr[i];
            answer++;
        }
        return answer;
    }
}