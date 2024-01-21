class Solution {
    public String solution(String my_string) {
        String answer = "";
        for(char c : my_string.toCharArray()){
            if('a' <= c && c <= 'z') answer += (char)(c - 'a' + 'A');
            else answer += (char)(c - 'A' + 'a');
            
        }
        return answer;
    }
}