import java.util.*;
class Solution {
    public String solution(String rsp) {
        HashMap<Character, Character> map = new HashMap<>();
        map.put('2', '0');
        map.put('0', '5');
        map.put('5', '2');
        String answer = "";
        for(char c : rsp.toCharArray()) answer += map.get(c);
        return answer;
    }
}