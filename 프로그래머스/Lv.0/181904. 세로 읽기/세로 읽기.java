class Solution {
    public String solution(String my_string, int m, int c) {
        String answer = "";
        int y = 0;
        int pos = c - 1;
        while(pos < my_string.length()){
            answer += my_string.charAt(pos);
            y++;
            pos = (y * m) + (c - 1);
        }
        return answer;
    }
}