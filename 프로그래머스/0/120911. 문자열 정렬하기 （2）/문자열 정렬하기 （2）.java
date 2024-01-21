import java.util.*;
import java.util.stream.Collectors;
class Solution {
    public String solution(String my_string) {
        char[] chars = my_string.toCharArray();
        for(int i = 0; i < chars.length; i++){
            if('A' <= chars[i] && chars[i] <= 'Z') chars[i] = (char)(chars[i] - 'A' + 'a');
        }
        Arrays.sort(chars);
        return String.valueOf(chars);
    }
}