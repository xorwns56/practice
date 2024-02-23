import java.util.*;
class Solution {
    public String solution(String[] participant, String[] completion) {
        HashMap<String, Integer> map = new HashMap<>();
        for(int j = 0; j < participant.length; j++){
            if(map.containsKey(participant[j])) map.put(participant[j], map.get(participant[j]) + 1);
            else map.put(participant[j], 1);
        }
        for(int j = 0; j < completion.length; j++){
            if(map.get(completion[j]) == 1) map.remove(completion[j]);
            else map.put(completion[j], map.get(completion[j]) - 1);
        }
        String answer = "";
        for(String key : map.keySet())  answer = key;
        return answer;
    }
}