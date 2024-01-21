class Solution {
    public int solution(String s) {
        String[] sp = s.split("\\s");
        int answer = 0;
        for(int i = 0; i < sp.length; i++){
            if(sp[i].equals("Z")) answer -= Integer.parseInt(sp[i - 1]);
            else answer += Integer.parseInt(sp[i]);;
        }
        return answer;
    }
}