import java.util.*;
class Solution {
    public int[] solution(String[] keymap, String[] targets) {
        HashMap<Character, Integer> map = new HashMap<>();
        for(int i = 0; i < keymap.length; i++){
            char[] chars = keymap[i].toCharArray();
            for(int j = 0; j < chars.length; j++){
                int touch = j + 1;
                if(map.containsKey(chars[j])) touch = Math.min(touch, map.get(chars[j]));
                map.put(chars[j], touch);
            }
        }
        int[] answer = new int[targets.length];
        for(int i = 0; i < targets.length; i++){
            int touches = 0;
            char[] chars = targets[i].toCharArray();
            for(int j = 0; j < chars.length; j++){
                if(!map.containsKey(chars[j])){
                    touches = -1;
                    break;
                }else touches += map.get(chars[j]);
            }
            answer[i] = touches;
        }
        return answer;
    }
}