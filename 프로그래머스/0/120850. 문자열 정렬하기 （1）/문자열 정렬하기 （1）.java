import java.util.*;
class Solution {
    public int[] solution(String my_string) {
        return my_string.chars().filter(c->'0'<=c&&c<='9').map(i->i-'0').sorted().toArray();
    }
}