import java.util.*;
class Solution {
    public String solution(String my_string) {
        HashSet<Character> set = new HashSet<>();
        String answer = "";
        for(char c : my_string.toCharArray()){
            if(!set.contains(c)){
                answer += c;
                set.add(c);
            }
        }
        return answer;
    }
}