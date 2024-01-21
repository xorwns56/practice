import java.util.*;
class Solution {
    public int solution(String before, String after) {
        char[] before_chars = before.toCharArray();
        char[] after_chars = after.toCharArray();
        Arrays.sort(before_chars);
        Arrays.sort(after_chars);
        return String.valueOf(before_chars).equals(String.valueOf(after_chars)) ? 1 : 0;
    }
}