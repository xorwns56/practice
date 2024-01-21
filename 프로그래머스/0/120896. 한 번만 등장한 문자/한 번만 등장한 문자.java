import java.util.*;
class Solution {
    public String solution(String s) {
        TreeMap<Character, Integer> map = new TreeMap<>();
        for(char c : s.toCharArray()) map.put(c, map.getOrDefault(c, 0) + 1);
        String answer = "";
        for(Map.Entry<Character, Integer> entry : map.entrySet()){
            if(entry.getValue() == 1) answer += entry.getKey();
        }
        return answer;
    }
}