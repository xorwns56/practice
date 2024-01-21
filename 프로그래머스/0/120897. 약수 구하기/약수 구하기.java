import java.util.*;
class Solution {
    public int[] solution(int n) {
        LinkedHashSet<Integer> set = new LinkedHashSet<>();
        for(int i = 1; i <= n; i++){
            if(n % i == 0) set.add(i);
        }
        int[] answer = new int[set.size()];
        Iterator<Integer> iter = set.iterator();
        int i = 0;
        while(iter.hasNext()) answer[i++] = iter.next();
        return answer;
    }
}