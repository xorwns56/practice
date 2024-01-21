import java.util.*;
class Solution {
    public int[] solution(int n) {
        List<Integer> list = new ArrayList<>();
        int rest = n;
        for(int i = 2; i <= n; i++){
            boolean prime = true;
            for(int j = 2; j <= (int)Math.sqrt(i); j++){
                if(i % j == 0){
                    prime = false;
                    break;
                }
            }
            if(rest % i == 0) list.add(i);
            while(rest % i == 0) rest /= i;
        }
        int[] answer = new int[list.size()];
        for(int i = 0; i < answer.length; i++) answer[i] = list.get(i);
        return answer;
    }
}