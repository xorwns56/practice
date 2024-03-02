import java.util.*;
class Solution {
    public void put(HashMap<String, Integer> map, char c1, char c2){
        if('a' <= c1 && c1 <= 'z') c1 += 'A' - 'a';
        else if(!('A' <= c1 && c1 <= 'Z')) return;
        if('a' <= c2 && c2 <= 'z') c2 += 'A' - 'a';
        else if(!('A' <= c2 && c2 <= 'Z')) return;
        String key = c1 + "" + c2;
        map.put(key, map.getOrDefault(key, 0) + 1);
    }
    public int solution(String str1, String str2) {
        char[] chars1 = str1.toCharArray();
        char[] chars2 = str2.toCharArray();
        HashMap<String, Integer> map1 = new HashMap<>();
        HashMap<String, Integer> map2 = new HashMap<>();
        for(int i = 1; i < chars1.length; i++) put(map1, chars1[i - 1], chars1[i]);
        for(int i = 1; i < chars2.length; i++) put(map2, chars2[i - 1], chars2[i]);
        int min_sum = 0;
        int max_sum = 0;
        for(String key : map1.keySet()){
            int val1 = map1.getOrDefault(key, 0);
            int val2 = 0;
            if(map2.containsKey(key)){
                val2 = map2.get(key);
                map2.remove(key);
            }
            min_sum += Math.min(val1, val2);
            max_sum += Math.max(val1, val2);
        }
        for(Integer val : map2.values()) max_sum += val;
        return (int)((min_sum != max_sum ? (double)min_sum / max_sum : 1) * 65536);
    }
}