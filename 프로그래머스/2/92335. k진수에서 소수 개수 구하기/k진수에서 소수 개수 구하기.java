import java.util.*;
class Solution {
    public int solution(int n, int k) {
        List<Long> list = new ArrayList<>();
        long tmp = 0;
        int count = 0;
        while(n > 0){
            if(n % k == 0){
                list.add(tmp);
                tmp = count = 0;
            }else tmp += (n % k) * Math.pow(10, count++);
            n /= k;
        }
        if(tmp > 0) list.add(tmp);
        int answer = 0;
        for(int i = 0; i < list.size(); i++){
            boolean isPrime = true;
            if(list.get(i) <= 1) isPrime = false;
            for(int j = 2; j <= (int)Math.sqrt(list.get(i)); j++){
                if(list.get(i) % j == 0){
                    isPrime = false;
                    break;
                }
            }
            if(isPrime) answer++;
        }
        return answer;
    }
    
}