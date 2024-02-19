import java.util.*;
class Solution {
    public int solution(int[] d, int budget) {
        Arrays.sort(d);
        int i;
        for(i = 0; i < d.length; i++){
            budget -= d[i];
            if(budget < 0) break;
        }
        return i;
    }
}