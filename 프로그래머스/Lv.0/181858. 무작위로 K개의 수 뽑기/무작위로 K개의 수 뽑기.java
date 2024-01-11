import java.util.*;
class Solution {
    public int[] solution(int[] arr, int k) {
        HashSet<Integer> hashSet = new HashSet<>();
        int[] answer = new int[k];
        int idx = 0;
        for(int i = 0; i < arr.length; i++){
            if(hashSet.add(arr[i]) && idx < k) answer[idx++] = arr[i];
        }
        for(int i = idx; i < answer.length; i++) answer[i] = -1;
        return answer;
    }
}