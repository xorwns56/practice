import java.util.*;
class Solution {
    public int solution(int[] topping) {
        HashMap<Integer, Integer> map1 = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> map2 = new HashMap<Integer, Integer>();
        for(int i = 0; i < topping.length; i++) map2.put(topping[i], map2.getOrDefault(topping[i], 0) + 1);
        int answer = 0;
        for(int i = 0; i < topping.length; i++){
            map1.put(topping[i], map1.getOrDefault(topping[i], 0) + 1);
            int map2_count = map2.get(topping[i]);
            if(map2_count == 1) map2.remove(topping[i]);
            else map2.put(topping[i], map2_count - 1);
            if(map1.size() == map2.size()) answer++;
        }
        return answer;
    }
}