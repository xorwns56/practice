import java.util.*;
class Solution {
    public int[] solution(int[] arr) {
        int min = 0;
        for(int i = 1; i < arr.length; i++){
            if(arr[i] < arr[min]) min = i;
        }
        int[] answer = new int[arr.length - 1];
        int idx = 0;
        for(int i = 0; i < arr.length; i++){
            if(i == min) continue;
            answer[idx++] = arr[i];
        }
        return answer.length == 0 ? new int[] { -1 } : answer;
    }
}