import java.util.*;
class Solution {
    public String[] solution(String[] str_list) {
        for(int i = 0; i < str_list.length; i++){
            if(str_list[i].charAt(0) == 'l') return Arrays.copyOfRange(str_list, 0, i);
            else if(str_list[i].charAt(0) == 'r') return Arrays.copyOfRange(str_list, i + 1, str_list.length);
        }
        return new String[]{};
    }
}