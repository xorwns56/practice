import java.util.*;
class Solution {
    public int solution(int[] elements) {
        HashSet<Integer> set = new HashSet<>();
        for(int i = 1; i <= elements.length; i++){
            int moving_sum = 0;
            for(int j = 0; j < i; j++) moving_sum += elements[j];
            for(int j = 0; j < elements.length; j++){
                set.add(moving_sum);
                moving_sum -= elements[j];
                moving_sum += elements[(j + i) % elements.length];
            }
        }
        return set.size();
    }
}