import java.util.*;
class Solution {
    public int[] solution(String msg) {
        HashMap<String, Integer> map = new HashMap<>();
        int key_index = 1;
        for(char c = 'A'; c <= 'Z'; c++) map.put(Character.toString(c), key_index++);
        char[] chars = msg.toCharArray();
        String str = "";
        List<Integer> list = new ArrayList<>();
        for(int i = 0; i < chars.length; i++){
            str += chars[i];
            if(map.get(str) == null){
                list.add(map.get(str.substring(0, str.length() - 1)));
                map.put(str, key_index++);
                str = "";
                i--;
            }
        }
        if(!str.isEmpty()) list.add(map.get(str));
        int[] answer = new int[list.size()];
        for(int i = 0; i < list.size(); i++) answer[i] = list.get(i);
        return answer;
    }
}