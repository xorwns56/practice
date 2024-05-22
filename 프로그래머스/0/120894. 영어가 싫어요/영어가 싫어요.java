import java.util.*;
class Solution {
    public long solution(String numbers) {
        HashMap<String, Integer> map = new HashMap<>();
        String[] numStr = { "zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine" };
        for(int i = 0; i < numStr.length; i++) map.put(numStr[i], i);
        String tmp = "";
        long answer = 0;
        for(int i = 0; i < numbers.length(); i++){
            tmp += numbers.charAt(i);
            if(map.containsKey(tmp)){
                answer = answer * 10 + map.get(tmp);
                tmp = "";
            }
        }
        return answer;
    }
}