import java.util.*;
class Solution {
    public int solution(String[] want, int[] number, String[] discount) {
        int answer = 0;
        HashMap<String, Integer> map = new HashMap<>();
        for(int i = 0; i < discount.length; i++){
            map.put(discount[i], map.getOrDefault(discount[i], 0) + 1);
            if(i >= 9){
                boolean discount_all = true;
                for(int j = 0; j < want.length; j++){
                    if(map.getOrDefault(want[j], 0) < number[j]) discount_all = false;
                }
                if(discount_all) answer++;
                map.put(discount[i - 9], map.getOrDefault(discount[i - 9], 0) - 1);
            }
        }
        return answer;
    }
}