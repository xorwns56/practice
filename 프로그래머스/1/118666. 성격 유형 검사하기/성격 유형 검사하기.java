import java.util.*;
class Solution {
    public String solution(String[] survey, int[] choices) {
        String answer = "";
        char[] type_char = new char[] { 'R', 'T', 'C', 'F', 'J', 'M', 'A', 'N'};
        int[] type_score = new int[type_char.length];
        HashMap<Character, Integer> map = new HashMap<>();
        for(int i = 0; i < type_char.length; i++) map.put(type_char[i], i);
        for(int i = 0; i < survey.length; i++){
            int score = choices[i] - 4;
            type_score[map.get(survey[i].charAt(score < 0 ? 0 : 1))] += Math.abs(score);
        }
        for(int i = 0; i < type_score.length; i += 2){
            answer += type_score[i] >= type_score[i + 1] ? type_char[i] : type_char[i + 1];
        }
        return answer;
    }
}