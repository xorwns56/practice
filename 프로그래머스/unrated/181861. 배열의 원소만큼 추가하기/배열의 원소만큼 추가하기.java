import java.util.*;
class Solution {
    public int[] solution(int[] arr) {
        LinkedList<Integer> list = new LinkedList<>();
        for(int i = 0; i < arr.length; i++){
             for(int j = 0; j < arr[i]; j++) list.add(arr[i]);
        }
        return list.stream().mapToInt(i -> i).toArray();
    }
}