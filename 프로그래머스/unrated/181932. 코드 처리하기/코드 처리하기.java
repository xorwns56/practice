class Solution {
    public String solution(String code) {
        int mode = 0;
        String answer = "";
        for(int i = 0; i < code.length(); i++){
            char c = code.charAt(i);
            if(c == '1') mode = mode ^ 1;
            else if((i & 1) == mode) answer += c;
        }
        return answer.isEmpty() ? "EMPTY" : answer;
    }
}