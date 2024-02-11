import java.util.*;
class Solution {
    public long solution(long n) {
        List<Long> list = new ArrayList<>();
        while(n > 0){
            list.add(n % 10);
            n /= 10;
        }
        Collections.sort(list);
        long answer = 0;
        for(int i = list.size() - 1; i >= 0; i--) answer = answer * 10 + list.get(i);
        return answer;
    }
}