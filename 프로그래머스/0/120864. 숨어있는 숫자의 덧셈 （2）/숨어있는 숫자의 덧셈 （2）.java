class Solution {
    public int solution(String my_string) {
        int answer = 0;
        int tmp = 0;
        for(char c : my_string.toCharArray()){
            if('0' <= c && c <= '9'){
                if(tmp > 0) tmp = tmp * 10;
                tmp += c - '0';
            }else{
                answer += tmp;
                tmp = 0;
            }
        }
        if(tmp > 0) answer += tmp;
        return answer;
    }
}