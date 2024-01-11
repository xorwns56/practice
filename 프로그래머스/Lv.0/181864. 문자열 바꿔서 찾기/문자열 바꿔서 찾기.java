import java.util.*;
class Solution {
    public int solution(String myString, String pat) {
        char[] c = myString.toCharArray();
        for(int i = 0; i < c.length; i++) c[i] = c[i] == 'A' ? 'B' : 'A';
        return String.valueOf(c).contains(pat) ? 1 : 0;
    }
}